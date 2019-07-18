/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "porousZone.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

// adjust negative resistance values to be multiplier of max value
void Foam::porousZone::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = max(0, cmptMax(resist.value()));

    if (maxCmpt < 0)
    {
        FatalErrorIn
        (
            "Foam::porousZone::porousZone::adjustNegativeResistance"
            "(dimensionedVector&)"
        )   << "negative resistances! " << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt=0; cmpt < vector::nComponents; ++cmpt)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZone::porousZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    cellZoneID_(mesh_.cellZones().findZoneID(name)),
    coordSys_(dict, mesh),
    porosity_(1),
    C0_(0),
    C1_(0),
    D_("D", dimensionSet(0, -2, 0, 0, 0), tensor::zero),
    F_("F", dimensionSet(0, -1, 0, 0, 0), tensor::zero),
    rhoPri_(0),
    Cpri_(0),
    Maux_(0),
    Qaux_(0),
    Caux_(0),
    Taux_(0),
    Qepsilon_(0),
    TauxRelax_(1),
    nCellsAuxInlet_(0),
    firstCell_(0),
    auxUnitVector_(vector::zero),
    nVerticalCells_(0)
{
    Info<< "Creating porous zone: " << name_ << endl;

    bool foundZone = (cellZoneID_ != -1);
    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const fvMesh&, const word&, const dictionary&)"
        )   << "cannot find porous cellZone " << name_
            << exit(FatalError);
    }


    // porosity
    if (dict_.readIfPresent("porosity", porosity_))
    {
        if (porosity_ <= 0.0 || porosity_ > 1.0)
        {
            FatalIOErrorIn
            (
                "Foam::porousZone::porousZone"
                "(const fvMesh&, const word&, const dictionary&)",
                dict_
            )
                << "out-of-range porosity value " << porosity_
                << exit(FatalIOError);
        }
    }

    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        dictPtr->readIfPresent("C0", C0_);
        dictPtr->readIfPresent("C1", C1_);
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        // local-to-global transformation tensor
        const tensor& E = coordSys_.R();

        dimensionedVector d(vector::zero);
        if (dictPtr->readIfPresent("d", d))
        {
            if (D_.dimensions() != d.dimensions())
            {
                FatalIOErrorIn
                (
                    "Foam::porousZone::porousZone"
                    "(const fvMesh&, const word&, const dictionary&)",
                    dict_
                )   << "incorrect dimensions for d: " << d.dimensions()
                    << " should be " << D_.dimensions()
                    << exit(FatalIOError);
            }

            adjustNegativeResistance(d);

            D_.value().xx() = d.value().x();
            D_.value().yy() = d.value().y();
            D_.value().zz() = d.value().z();
            D_.value() = (E & D_ & E.T()).value();
        }

        dimensionedVector f(vector::zero);
        if (dictPtr->readIfPresent("f", f))
        {
            if (F_.dimensions() != f.dimensions())
            {
                FatalIOErrorIn
                (
                    "Foam::porousZone::porousZone"
                    "(const fvMesh&, const word&, const dictionary&)",
                    dict_
                )   << "incorrect dimensions for f: " << f.dimensions()
                    << " should be " << F_.dimensions()
                    << exit(FatalIOError);
            }

            adjustNegativeResistance(f);

            // leading 0.5 is from 1/2 * rho
            F_.value().xx() = 0.5*f.value().x();
            F_.value().yy() = 0.5*f.value().y();
            F_.value().zz() = 0.5*f.value().z();
            F_.value() = (E & F_ & E.T()).value();
        }
    }

    // provide some feedback for the user
    // writeDict(Info, false);

    // It is an error not to define anything
    if
    (
        C0_ <= VSMALL
     && magSqr(D_.value()) <= VSMALL
     && magSqr(F_.value()) <= VSMALL
    )
    {
        FatalIOErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const fvMesh&, const word&, const dictionary&)",
            dict_
        )   << "neither powerLaw (C0/C1) "
               "nor Darcy-Forchheimer law (d/f) specified"
            << exit(FatalIOError);
    }

    // Heat rate value and temperature of heat exchanger
    if (const dictionary* dictPtr = dict_.subDictPtr("heatTransfer"))
    {
        Info<< "Reading porous heatTransfer: Maux and Taux" << nl;
        dictPtr->lookup("rhoPri") >> rhoPri_;
        dictPtr->lookup("Cpri") >> Cpri_;

        dictPtr->lookup("Maux") >> Maux_;
        dictPtr->lookup("Qaux") >> Qaux_;
        dictPtr->lookup("Caux") >> Caux_;
        dictPtr->lookup("Taux") >> Taux_;
        dictPtr->lookup("Qepsilon") >> Qepsilon_;

        dictPtr->lookup("firstCell") >> firstCell_;
        dictPtr->lookup("auxUnitVector") >> auxUnitVector_;
        dictPtr->lookup("nVerticalCells") >> nVerticalCells_;
        dictPtr->lookup("nCellsAuxInlet") >> nCellsAuxInlet_;

        dictPtr->lookup("TauxRelax") >> TauxRelax_;

        Info<< "Maux = " << Maux_
            << ", Taux = " << Taux_
            << ", Caux = " << Caux_
            << ", nCellsAuxInlet = " << nCellsAuxInlet_
            << ", Qepsilon = " << Qepsilon_
            << endl;

        if (mag(auxUnitVector_) < SMALL)
        {
            FatalIOErrorInFunction(dict)
                << "auxUnitVector for heat transfer zone "
                << name_ << " has zero length"
                << exit(FatalIOError);
        }
        else
        {
            auxUnitVector_ /= mag(auxUnitVector_);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZone::addResistance(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                Udiag,
                cells,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                Udiag,
                cells,
                V,
                geometricOneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }
}


void Foam::porousZone::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;

    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                AU,
                cells,
                geometricOneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is
        // interpolated for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


void Foam::porousZone::addHeatResistance
(
    fvScalarMatrix& hTEqn,
    const volScalarField& T,
    volScalarField& Taux,
    volScalarField& Qaux,
    const volVectorField& U,
    const volScalarField& Macro,
    const volScalarField& posFlux
) const
{
    if (cellZoneID_ == -1 && mag(Maux_) < SMALL)
    {
        return;
    }

    scalarField& QSource = hTEqn.source();
    const scalarField& Ti = T.internalField();
    scalarField& Tauxi = Taux.internalField();
    scalarField& Qauxi = Qaux.internalField();
    const vectorField& Ui = U.internalField();
    const scalarField& Macroi = Macro.internalField();
    const scalarField& posFluxi = posFlux.internalField();



    if (hTEqn.dimensions() == dimensionSet(1, -2, -3, 0, 0))
    {
        addHeatSource
        (
            Macroi,
            posFluxi,
            QSource,
            Qauxi,
            Ti,
            Tauxi,
            Ui,
            mesh_.lookupObject<volScalarField>("rho")
        );
    }
    else if (hTEqn.dimensions() == dimensionSet(0, 3, -1, 1, 0))
    {
        addHeatSource
        (
            Macroi,
            posFluxi,
            QSource,
            Qauxi,
            Ti,
            Tauxi,
            Ui,
            geometricOneField()
        );
    }
    else
    {
        Info<< "No hEqn or TEqn, exiting" << nl;
        return;
    }
}

void Foam::porousZone::macroCellOrder
(
    volScalarField& Taux,
    volScalarField& Macro,
    volScalarField& posFlux,
    const surfaceScalarField& phi
) const
{
    Info<< "Creating cellsOrdered list for porous heat transfer zone "
        << name_ << nl << endl;

    // Get current zone
    const cellZone& curZone = mesh_.cellZones()[cellZoneID_];

    // Get cell centres
    const vectorField& cellCentres = mesh_.C().internalField();

    // Get cell-cell addressing
    const labelListList& cc = mesh_.cellCells();

    scalarField& Tauxi = Taux.internalField();

    // Get inlet cells
    labelList cellsAuxInlet(nCellsAuxInlet_, -1);

    if (nCellsAuxInlet_ > 0)
    {
        cellsAuxInlet[0] = firstCell_;
    }
    else
    {
        WarningInFunction
            << "Number of aux inlet cells is set to zero: " << nCellsAuxInlet_
            << "Reconsider the definition of macro cell order"
            << endl;
    }

    // Creating horizontal list of cells, cellsAuxInlet[nCellsAuxInlet_]
    label newCounter = 1;

    // Note
    // Terrible search algorithms
    // Will work only in serial.  Rewrite required
    // HJ, 17/Jul/2019

    forAll (cellsAuxInlet, i)
    {
        const labelList& cellNb = cc[cellsAuxInlet[i]];

        forAll (cellNb, ii)
        {
            bool isInsideNb = false;

            forAll (cellsAuxInlet, iii)
            {
                if (cellNb[ii] == cellsAuxInlet[iii])
                {
                    isInsideNb = true;
                    break;
                }
            }

            if (!isInsideNb)
            {
                if (curZone.whichCell(cellNb[ii]) != -1)
                {
                    if
                    (
                        mag
                        (
                            (
                                (
                                    cellCentres[cellsAuxInlet[i]]
                                  - cellCentres[cellNb[ii]]
                                )/
                                mag
                                (
                                    cellCentres[cellsAuxInlet[i]]
                                  - cellCentres[cellNb[ii]]
                                )
                            )
                          & auxUnitVector_
                        ) < 0.5
                    )
                    {
                        cellsAuxInlet[newCounter] = cellNb[ii];
                        ++newCounter;
                    }
                }
            }
        }
    }

    scalarField& Macroi = Macro.internalField();

    labelList cellsOrdered(curZone.size(), -1);

    // Writing horizontal cells into list cellsOrdered[cells.size()]
    forAll (cellsAuxInlet, i)
    {
        cellsOrdered[i*nVerticalCells_] = cellsAuxInlet[i];
    }

    // Creating cellsOrdered list of ordered horizontal cells
    label counter = 1;

    forAll (cellsOrdered, i)
    {
        Macroi[cellsOrdered[i]] = i;
        Tauxi[cellsOrdered[i]] = Taux_;

        if ((i > 1) && (i % nVerticalCells_ == 0))
        {
            ++counter;
        }

        const labelList& cellNb =cc[cellsOrdered[i]];

        forAll (cellNb, iii)
        {
            bool isInsideNb = false;

            forAll (cellsOrdered, iiii)
            {
                if (cellNb[iii] == cellsOrdered[iiii])
                {
                    isInsideNb = true;
                    break;
                }
            }

            if (!isInsideNb)
            {
                if (curZone.whichCell(cellNb[iii]) != -1)
                {
                    if
                    (
                        mag
                        (
                            (
                                (
                                    cellCentres[cellsOrdered[i]]
                                  - cellCentres[cellNb[iii]]
                                )/
                                mag
                                (
                                    cellCentres[cellsOrdered[i]]
                                  - cellCentres[cellNb[iii]]
                                )
                            )
                          & auxUnitVector_
                        ) > 0.85
                    )
                    {
                        cellsOrdered[counter] = cellNb[iii];
                        ++counter;
                    }
                }
            }
        }
    }

    scalarField& posFluxi = posFlux.internalField();

    // Calculating mass flow through each cell
    forAll (cellsOrdered, i)
    {
        const labelList& cellFaces = mesh_.cells()[cellsOrdered[i]];

        forAll (cellFaces, ii)
        {
            label faceI = cellFaces[ii];

            if (mesh_.isInternalFace(faceI))
            {
                if (mesh_.faceOwner()[faceI] == cellsOrdered[i])
                {
                    if (phi.internalField()[faceI] > 0.0)
                    {
                        posFluxi[cellsOrdered[i]] += phi.internalField()[faceI];
                    }
                }
                else
                {
                    if (phi.internalField()[faceI] < 0.0)
                    {
                        posFluxi[cellsOrdered[i]] +=
                            -phi.internalField()[faceI];
                    }
                }
            }
            else
            {
                const label patchI = mesh_.boundaryMesh().whichPatch(faceI);

                const label faceIL =
                    mesh_.boundaryMesh()[patchI].whichFace(faceI);

                if (phi.boundaryField()[patchI][faceIL] > 0.0)
                {
                    posFluxi[cellsOrdered[i]] +=
                        phi.boundaryField()[patchI][faceIL];
                }
            }
        }
    }
}


void Foam::porousZone::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("name") << zoneName() << token::END_STATEMENT << nl;
    }
    else
    {
        os  << indent << zoneName() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    if (dict_.found("note"))
    {
        os.writeKeyword("note") << string(dict_.lookup("note"))
            << token::END_STATEMENT << nl;
    }

    coordSys_.writeDict(os, true);

    if (dict_.found("porosity"))
    {
        os.writeKeyword("porosity")
            << porosity() << token::END_STATEMENT
            << nl;
    }

    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        os << indent << "powerLaw";
        dictPtr->write(os);
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        os << indent << "Darcy";
        dictPtr->write(os);
    }

    // Heat transfer inputs
    if (const dictionary* dictPtr = dict_.subDictPtr("heatTransfer"))
    {
        os << indent << "heatTransfer";
        dictPtr->write(os);
    }

    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const porousZone& pZone)
{
    pZone.writeDict(os);
    return os;
}

// ************************************************************************* //
