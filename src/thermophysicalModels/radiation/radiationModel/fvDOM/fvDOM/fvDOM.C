/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "fvDOM.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "mathematicalConstants.H"
#include "radiationConstants.H"

using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(fvDOM, 0);
        addToRadiationRunTimeSelectionTables(fvDOM);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::fvDOM::initialise()
{
    if (mesh().nSolutionD() == 3)    //3D
    {
        nRay_ = 4*nPhi_*nTheta_;
        IRay_.setSize(nRay_);
        scalar deltaPhi = pi/(2.0*nPhi_);
        scalar deltaTheta = pi/nTheta_;
        label i = 0;
        for (label n = 1; n <= nTheta_; n++)
        {
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar thetai = (2.0*n - 1.0)*deltaTheta/2.0;
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh(),
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }
    else
    {
        if (mesh().nSolutionD() == 2)    //2D (X & Y)
        {
            scalar thetai = piByTwo;
            scalar deltaTheta = pi;
            nRay_ = 4*nPhi_;
            IRay_.setSize(nRay_);
            scalar deltaPhi = pi/(2.0*nPhi_);
            label i = 0;
            for (label m = 1; m <= 4*nPhi_; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh(),
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
        else    //1D (X)
        {
            scalar thetai = piByTwo;
            scalar deltaTheta = pi;
            nRay_ = 2;
            IRay_.setSize(nRay_);
            scalar deltaPhi = pi;
            label i = 0;
            for (label m = 1; m <= 2; m++)
            {
                scalar phii = (2.0*m - 1.0)*deltaPhi/2.0;
                IRay_.set
                (
                    i,
                    new radiativeIntensityRay
                    (
                        *this,
                        mesh(),
                        phii,
                        thetai,
                        deltaPhi,
                        deltaTheta,
                        nLambda_,
                        absorptionEmission_,
                        blackBody_,
                        i
                    )
                );
                i++;
            }
        }
    }


    // Construct absorption field for each wavelength
    forAll(aLambda_, lambdaI)
    {
        aLambda_.set
        (
            lambdaI,
            new volScalarField
            (
                IOobject
                (
                    "aLambda_" + Foam::name(lambdaI) ,
                    mesh().time().timeName(),
                    T_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                a_
            )
        );

        Qin_.set
        (
            lambdaI,
            volScalarField::GeometricBoundaryField
            (
                mesh().boundary(),
                mesh().V(),           // Dummy internal field,
                calculatedFvPatchScalarField::typeName
            )
        );
        Qin_[lambdaI] = 0;

        Qem_.set
        (
            lambdaI,
            volScalarField::GeometricBoundaryField
            (
                mesh().boundary(),
                mesh().V(),           // Dummy internal field,
                calculatedFvPatchScalarField::typeName
            )
        );
        Qem_[lambdaI] = 0;
    }

    Info<< "fvDOM : Allocated " << IRay_.size()
        << " rays with average orientation:" << nl;
    forAll (IRay_, i)
    {
        Info<< '\t' << IRay_[i].I().name()
            << '\t' << IRay_[i].dAve() << nl;
    }


    if (cacheDiv_)
    {
        Info<< "Caching div fvMatrix..."<< endl;
        for (label lambdaI = 0; lambdaI < nLambda_; lambdaI++)
        {
            fvRayDiv_[lambdaI].setSize(nRay_);

            forAll(IRay_, rayId)
            {
                const surfaceScalarField Ji(IRay_[rayId].dAve() & mesh().Sf());
                volScalarField& iRayLambdaI = IRay_[rayId].ILambda(lambdaI);

                fvRayDiv_[lambdaI].set
                (
                    rayId,
                    new fvScalarMatrix
                    (
                        fvm::div(Ji, iRayLambdaI, "div(Ji,Ii_h)")
                    )
                );
            }
        }
    }

    forAll(IRay_, rayId)
    {
        if (omegaMax_ <  IRay_[rayId].omega())
        {
            omegaMax_ = IRay_[rayId].omega();
        }
        Info<< '\t' << IRay_[rayId].I().name() << " : " << "omega : "
            << '\t' << IRay_[rayId].omega() << nl;
    }

    Info<< endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::fvDOM(const volScalarField& T)
:
    radiationModel(typeName, T),
    G_
    (
        IOobject
        (
            "G",
            mesh().time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            T.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh().time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    Qem_(nLambda_),
    Qin_(nLambda_),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    fvRayDiv_(nLambda_),
    cacheDiv_(coeffs_.lookupOrDefault<bool>("cacheDiv", false)),
    omegaMax_(0)
{
    initialise();
}


Foam::radiation::fvDOM::fvDOM(const word& type, const volScalarField& T)
:
    radiationModel(type, T),
    G_
    (
        IOobject
        (
            "G",
            mesh().time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            T.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh().time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    Qem_(nLambda_),
    Qin_(nLambda_),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    fvRayDiv_(nLambda_),
    cacheDiv_(coeffs_.lookupOrDefault<bool>("cacheDiv", false)),
    omegaMax_(0)
{
    initialise();
}


Foam::radiation::fvDOM::fvDOM
(
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(typeName, dict, T),
    G_
    (
        IOobject
        (
            "G",
            mesh().time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("G", dimMass/pow3(dimTime), 0.0)
    ),
    Qr_
    (
        IOobject
        (
            "Qr",
            mesh_.time().timeName(),
            T.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
    ),
    a_
    (
        IOobject
        (
            "a",
            mesh().time().timeName(),
            T.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("a", dimless/dimLength, 0.0)
    ),
    nTheta_(readLabel(coeffs_.lookup("nTheta"))),
    nPhi_(readLabel(coeffs_.lookup("nPhi"))),
    nRay_(0),
    nLambda_(absorptionEmission_->nBands()),
    aLambda_(nLambda_),
    blackBody_(nLambda_, T),
    IRay_(0),
    Qem_(nLambda_),
    Qin_(nLambda_),
    convergence_(coeffs_.lookupOrDefault<scalar>("convergence", 0.0)),
    maxIter_(coeffs_.lookupOrDefault<label>("maxIter", 50)),
    fvRayDiv_(nLambda_),
    cacheDiv_(coeffs_.lookupOrDefault<bool>("cacheDiv", false)),
    omegaMax_(0)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fvDOM::~fvDOM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::GeometricBoundaryField>
Foam::radiation::fvDOM::Qin() const
{
    tmp<volScalarField::GeometricBoundaryField> tQin
    (
        new volScalarField::GeometricBoundaryField
        (
            mesh().boundary(),
            mesh().V(),           // Dummy internal field,
            calculatedFvPatchScalarField::typeName
        )
    );
    volScalarField::GeometricBoundaryField& sumQin = tQin();

    sumQin = 0;

    forAll(Qin_, lambdaI)
    {
        sumQin += Qin(lambdaI);
    }

    return tQin;
}


Foam::tmp<Foam::volScalarField::GeometricBoundaryField>
Foam::radiation::fvDOM::Qem() const
{
    tmp<volScalarField::GeometricBoundaryField> tsumQem
    (
        new volScalarField::GeometricBoundaryField
        (
            mesh().boundary(),
            mesh().V(),           // Dummy internal field,
            calculatedFvPatchScalarField::typeName
        )
    );
    volScalarField::GeometricBoundaryField& sumQem = tsumQem();

    sumQem = 0;

    forAll(Qem_, lambdaI)
    {
        sumQem += Qem(lambdaI);
    }

    return tsumQem;
}


bool Foam::radiation::fvDOM::read()
{
    if (radiationModel::read())
    {
        // Only reading solution parameters - not changing ray geometry
        coeffs_.readIfPresent("convergence", convergence_);
        coeffs_.readIfPresent("maxIter", maxIter_);

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::fvDOM::calculate()
{
    absorptionEmission_->correct(a_, aLambda_);

    updateBlackBodyEmission();

    scalar maxResidual = 0;
    label radIter = 0;
    do
    {
        Info << "Radiation solver iter: " << radIter << endl;

        Info << "Updating Radiation BCs..." << flush;
        forAll(IRay_, rayI)
        {
            IRay_[rayI].updateBCs();
        }
        Info << "done." << endl;

        // For debug purposes, recalculate Qin and Qem to see if balances match
        if (debug)
        {
            // Update radiation balances
            forAll(aLambda_, lambdaI)
            {
                Qem_[lambdaI] = 0;
                Qin_[lambdaI] = 0;

                forAll(Qem_[lambdaI], patchI)
                {
                    // Loop over all rays
                    forAll(IRay_, rayI)
                    {
                        const fvPatchScalarField& curPatch =
                            IRay_[rayI].ILambda
                            (
                                lambdaI
                            ).boundaryField()[patchI];

                        scalarField incomingAngle =
                            curPatch.patch().nf() & IRay_[rayI].dAve();

                        Qin_[lambdaI][patchI] +=
                            pos(incomingAngle)*curPatch*incomingAngle;

                        Qem_[lambdaI][patchI] +=
                            neg(incomingAngle)*curPatch*(-incomingAngle);
                    }

                    Info << "Patch " << Qem_[lambdaI][patchI].patch().name()
                        << " band " << lambdaI
                        << ": Radiation incoming "
                        << sum(Qin_[lambdaI][patchI])
                        << "  outgoing "
                        << sum(Qem_[lambdaI][patchI])
                        << endl;
                }
            }
        }

        // Solve ray transport equations
        forAll(IRay_, rayI)
        {
            maxResidual = 0;
            scalar maxRayResidual = IRay_[rayI].correct();
            maxResidual = max(maxRayResidual, maxResidual);
        }

        // Update radiation balances
        forAll(aLambda_, lambdaI)
        {
            Qem_[lambdaI] = 0;
            Qin_[lambdaI] = 0;

            forAll(Qem_[lambdaI], patchI)
            {
                // Loop over all rays
                forAll(IRay_, rayI)
                {
                    const fvPatchScalarField& curPatch =
                        IRay_[rayI].ILambda(lambdaI).boundaryField()[patchI];

                    scalarField incomingAngle =
                        curPatch.patch().nf() & IRay_[rayI].dAve();

                    Qin_[lambdaI][patchI] +=
                        pos(incomingAngle)*curPatch*incomingAngle;

                    Qem_[lambdaI][patchI] +=
                        neg(incomingAngle)*curPatch*(-incomingAngle);
                }
            }
        }

        radIter++;
        Info << "Max residual: " << maxResidual << endl;

    } while(maxResidual > convergence_ && radIter < maxIter_);

    updateG();
}


Foam::tmp<Foam::volScalarField> Foam::radiation::fvDOM::Rp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Rp",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            4.0*a_*radiation::sigmaSB //absorptionEmission_->a()
        )
    );
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::radiation::fvDOM::Ru() const
{

    const DimensionedField<scalar, volMesh>& G =
        G_.dimensionedInternalField();
    const DimensionedField<scalar, volMesh> E =
        absorptionEmission_->ECont()().dimensionedInternalField();
    const DimensionedField<scalar, volMesh> a =
        a_.dimensionedInternalField();

    return  a*G - E;
}


void Foam::radiation::fvDOM::updateBlackBodyEmission()
{
    for (label j=0; j < nLambda_; j++)
    {
        blackBody_.correct(j, absorptionEmission_->bands(j));
    }
}


void Foam::radiation::fvDOM::updateG()
{
    G_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);
    Qr_ = dimensionedScalar("zero",dimMass/pow3(dimTime), 0.0);

    forAll(IRay_, rayI)
    {
        IRay_[rayI].addIntensity();
        G_ += IRay_[rayI].I()*IRay_[rayI].omega();
    }

    forAll(Qr_.boundaryField(), patchI)
    {
        forAll(aLambda_, lambdaI)
        {
            Qr_.boundaryField()[patchI] += Qin_[lambdaI][patchI] - Qem_[lambdaI][patchI];
        }
    }
}


void Foam::radiation::fvDOM::setRayIdLambdaId
(
    const word& name,
    label& rayId,
    label& lambdaId
) const
{
    // assuming name is in the form: CHARS_rayId_lambdaId
    size_type i1 = name.find_first_of("_");
    size_type i2 = name.find_last_of("_");

    rayId = readLabel(IStringStream(name.substr(i1+1, i2-1))());
    lambdaId = readLabel(IStringStream(name.substr(i2+1, name.size()-1))());
}


// ************************************************************************* //
