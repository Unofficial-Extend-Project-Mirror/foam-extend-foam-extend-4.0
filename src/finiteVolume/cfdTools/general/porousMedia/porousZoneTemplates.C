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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::porousZone::modifyDdt(fvMatrix<Type>& m) const
{
    if (porosity_ < 1)
    {
        const labelList& zoneCells = mesh_.cellZones()[cellZoneID_];

        forAll (zoneCells, i)
        {
            m.diag()[zoneCells[i]] *= porosity_;
            m.source()[zoneCells[i]] *= porosity_;
        }
    }
}


template<class RhoFieldType>
void Foam::porousZone::addPowerLawResistance
(
    scalarField& Udiag,
    const labelList& cells,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll (cells, i)
    {
        Udiag[cells[i]] +=
            V[cells[i]]*rho[cells[i]]*C0*pow(magSqr(U[cells[i]]), C1m1b2);
    }
}


template<class RhoFieldType>
void Foam::porousZone::addViscousInertialResistance
(
    scalarField& Udiag,
    vectorField& Usource,
    const labelList& cells,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    const tensor& D = D_.value();
    const tensor& F = F_.value();

    forAll (cells, i)
    {
        tensor dragCoeff = mu[cells[i]]*D + (rho[cells[i]]*mag(U[cells[i]]))*F;
        scalar isoDragCoeff = tr(dragCoeff);

        Udiag[cells[i]] += V[cells[i]]*isoDragCoeff;
        Usource[cells[i]] -=
            V[cells[i]]*((dragCoeff - I*isoDragCoeff) & U[cells[i]]);
    }
}


template<class RhoFieldType>
void Foam::porousZone::addPowerLawResistance
(
    tensorField& AU,
    const labelList& cells,
    const RhoFieldType& rho,
    const vectorField& U
) const
{
    const scalar C0 = C0_;
    const scalar C1m1b2 = (C1_ - 1.0)/2.0;

    forAll (cells, i)
    {
        AU[cells[i]] = AU[cells[i]]
          + I*(rho[cells[i]]*C0*pow(magSqr(U[cells[i]]), C1m1b2));
    }
}


template<class RhoFieldType>
void Foam::porousZone::addViscousInertialResistance
(
    tensorField& AU,
    const labelList& cells,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    const tensor& D = D_.value();
    const tensor& F = F_.value();

    forAll (cells, i)
    {
        AU[cells[i]] += mu[cells[i]]*D + (rho[cells[i]]*mag(U[cells[i]]))*F;
    }
}


template<class RhoFieldType>
void Foam::porousZone::addHeatSource
(
    const scalarField& Macro,
    const scalarField& posFlux,
    scalarField& QSource,
    scalarField& Qauxi,
    const scalarField& T,
    scalarField& Tauxi,
    const vectorField& U,
    const RhoFieldType& rho
) const
{
    // - Dual stream model (Single Effectiveness Model)
    const scalarField Tauxi_old = Tauxi; // create storePrevIter()
    const scalarField dT = Tauxi - T;

    const label nMacro(nCellsAuxInlet_*nVerticalCells_);
    scalarField Tmacro(nMacro, scalar(0));
    scalarField dTaux(nMacro, scalar(0));
    scalar QSum = 0;

    // Mass flow rate for individual channels
    const scalar qmAuxi = Maux_/nCellsAuxInlet_;

    // Get zone cells
    const labelList& zoneCells = mesh_.cellZones()[cellZoneID_];

    forAll (zoneCells, i)
    {
        scalar Qcell =
            Qepsilon_*rhoPri_*Cpri_*posFlux[zoneCells[i]]*dT[zoneCells[i]];

        // Heat in each macro (cell)
        Qauxi[zoneCells[i]] = Qcell;

        // make an int out of a macro
        const int macro = Macro[zoneCells[i]];

        // deltaTaux in each macro(cell)
        dTaux[macro] = Qcell/(Caux_*qmAuxi);

        // adding Heat to equation
        QSource[zoneCells[i]] += Qcell/(rhoPri_*Cpri_);

        // summing for total heat of HX
        QSum += Qcell;
    }

    reduce(dTaux, sumOp<scalarField>());

    Tmacro[0] = Taux_;

    forAll (Tmacro, i)
    {
        if (i > 0)
        {
            if (i % nVerticalCells_ == 0)
            {
                Tmacro[i] = Taux_;
            }
            else
            {
                Tmacro[i] = Tmacro[i - 1] - dTaux[i - 1];
            }
        }
    }

    reduce(QSum, sumOp<scalar>());

    // Adjust Taux_ to match the specified transferred heat

    scalar deltaTAux = (Qaux_ - QSum)/(Maux_*Caux_);

    Info<< "Heat exchanger: " << name_ << nl
        << "Q = " << QSum << nl
        << "Taux = " << Taux_ << nl
        << "delta Taux = " << deltaTAux
        << endl;

    // Relax the heat source
    forAll (zoneCells, i)
    {
        const int macro = Macro[zoneCells[i]];

        Tauxi[zoneCells[i]] = Tauxi_old[zoneCells[i]]
          + TauxRelax_*(Tmacro[macro] - Tauxi_old[zoneCells[i]]);
    }

    // Note: need better relaxation to speed up convergence close
    // to the matching value, when deltaTAux -> 0
    // HJ, 17/Jul/2019
    
    // Limit deltaTAux to 10% of Taux_
    deltaTAux = sign(deltaTAux)*Foam::min(0.1*Taux_, mag(deltaTAux));
    
    // Update Taux for next iteration
    Taux_ += TauxRelax_*deltaTAux;
}


// ************************************************************************* //
