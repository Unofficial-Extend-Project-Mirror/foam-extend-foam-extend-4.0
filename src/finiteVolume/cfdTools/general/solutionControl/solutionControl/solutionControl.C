/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "solutionControl.H"
#include "lduMatrix.H"
#include "demandDrivenData.H"
#include "fvc.H"
#include "slipFvPatchFields.H"
#include "symmetryFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solutionControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::solutionControl::read(const bool absTolOnly)
{
    const dictionary& solutionDict = this->dict();

    // Read solution controls
    nNonOrthCorr_ =
        solutionDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);
    momentumPredictor_ =
        solutionDict.lookupOrDefault("momentumPredictor", true);
    transonic_ = solutionDict.lookupOrDefault("transonic", false);
    consistent_ = solutionDict.lookupOrDefault("consistent", false);

    // Read residual information
    const dictionary residualDict
    (
        solutionDict.subOrEmptyDict("residualControl")
    );

    DynamicList<fieldData> data(residualControl_);

    forAllConstIter(dictionary, residualDict, iter)
    {
        const word& fName = iter().keyword();
        const label fieldI = applyToField(fName, false);
        if (fieldI == -1)
        {
            fieldData fd;
            fd.name = fName.c_str();

            if (absTolOnly)
            {
                fd.absTol = readScalar(residualDict.lookup(fName));
                fd.relTol = -1;
                fd.initialResidual = -1;
            }
            else
            {
                if (iter().isDict())
                {
                    const dictionary& fieldDict(iter().dict());
                    fd.absTol = readScalar(fieldDict.lookup("tolerance"));
                    fd.relTol = readScalar(fieldDict.lookup("relTol"));
                    fd.initialResidual = 0.0;
                }
                else
                {
                    FatalErrorIn
                    (
                        "void Foam::solutionControl::read"
                        "(const bool absTolOnly)"
                    )   << "Residual data for " << iter().keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }
            }

            data.append(fd);
        }
        else
        {
            fieldData& fd = data[fieldI];
            if (absTolOnly)
            {
                fd.absTol = readScalar(residualDict.lookup(fName));
            }
            else
            {
                if (iter().isDict())
                {
                    const dictionary& fieldDict(iter().dict());
                    fd.absTol = readScalar(fieldDict.lookup("tolerance"));
                    fd.relTol = readScalar(fieldDict.lookup("relTol"));
                }
                else
                {
                    FatalErrorIn
                    (
                        "void Foam::solutionControl::read"
                        "(const bool absTolOnly)"
                    )   << "Residual data for " << iter().keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }
            }
        }
    }

    residualControl_.transfer(data);

    if (debug)
    {
        forAll(residualControl_, i)
        {
            const fieldData& fd = residualControl_[i];
            Info<< "residualControl[" << i << "]:" << nl
                << "    name     : " << fd.name << nl
                << "    absTol   : " << fd.absTol << nl
                << "    relTol   : " << fd.relTol << nl
                << "    iniResid : " << fd.initialResidual << endl;
        }
    }
}


void Foam::solutionControl::read()
{
    read(false);
}


Foam::label Foam::solutionControl::applyToField
(
    const word& fieldName,
    const bool useRegEx
) const
{
    forAll(residualControl_, i)
    {
        if (useRegEx && residualControl_[i].name.match(fieldName))
        {
            return i;
        }
        else if (residualControl_[i].name == fieldName)
        {
            return i;
        }
    }

    return -1;
}


void Foam::solutionControl::storePrevIterFields() const
{
//    storePrevIter<label>();
    storePrevIter<scalar>();
    storePrevIter<vector>();
    storePrevIter<sphericalTensor>();
    storePrevIter<symmTensor>();
    storePrevIter<tensor>();
}


template<class Type>
void Foam::solutionControl::maxTypeResidual
(
    const word& fieldName,
    ITstream& data,
    scalar& firstRes,
    scalar& lastRes
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (mesh_.foundObject<fieldType>(fieldName))
    {
        const List<lduSolverPerformance> sp(data);
        firstRes = cmptMax(sp.first().initialResidual());
        lastRes = cmptMax(sp.last().initialResidual());
    }
}


Foam::scalar Foam::solutionControl::maxResidual
(
    const word& fieldName,
    ITstream& data,
    scalar& lastRes
) const
{
    scalar firstRes = 0;

    maxTypeResidual<scalar>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<vector>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<sphericalTensor>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<symmTensor>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<tensor>(fieldName, data, firstRes, lastRes);

    return firstRes;
}


const Foam::dimensionedScalar Foam::solutionControl::relaxFactor
(
    const Foam::volVectorField& U
) const
{
    scalar urf = 1;

    if (mesh_.solutionDict().relaxEquation(U.name()))
    {
        urf = mesh_.solutionDict().equationRelaxationFactor(U.name());
    }

    return dimensionedScalar("alphaU", dimless, urf);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionControl::solutionControl(fvMesh& mesh, const word& algorithmName)
:
    IOobject
    (
        "solutionControl",
        mesh.time().timeName(),
        mesh
    ),
    mesh_(mesh),
    residualControl_(),
    algorithmName_(algorithmName),
    nNonOrthCorr_(0),
    momentumPredictor_(true),
    transonic_(false),
    consistent_(false),
    corr_(0),
    corrNonOrtho_(0),
    aCoeffPtr_(NULL),
    faceUPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionControl::~solutionControl()
{
    deleteDemandDrivenData(aCoeffPtr_);
    deleteDemandDrivenData(faceUPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solutionControl::calcTimeConsistentFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& rAU,
    const fvVectorMatrix& ddtUEqn
) const
{
    // Store necessary fields for time consistency if they are not present
    if (!aCoeffPtr_ && !faceUPtr_)
    {
        aCoeffPtr_ = new volScalarField
        (
            IOobject
            (
                "aCoeff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0.0)
        );

        faceUPtr_ = new surfaceVectorField
        (
            IOobject
            (
                "faceU",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimVelocity, vector::zero)
        );
    }
    else if (!aCoeffPtr_ || !faceUPtr_)
    {
        FatalErrorIn
        (
            "tmp<surfaceScalarField> solutionControl::timeConsistentFlux"
            "\n("
            "\n    const volVectorField& U,"
            "\n    const volScalarField& rAU"
            "\n)"
        )   << "Either aCoeffPtr_ or faceUPtr_ is allocated while the"
            << " other is not. This must not happen."
            << endl
            << exit(FatalError);
    }

    // Now that we are sure that pointers are valid and the fields are there, we
    // can do the calculation of the flux, while storing necessary data for
    // future velocity reconstruction.
    // Algorithm:
    // 1. Update flux and aCoeff due to ddt discretisation
    // 2. Update flux and aCoeff due to under-relaxation
    // 3. Scale the flux with aCoeff, making sure that flux at fixed boundaries
    //    remains consistent

    // Get fields that will be updated
    volScalarField& aCoeff = *aCoeffPtr_;
    surfaceVectorField& faceU = *faceUPtr_;

    // Update face interpolated velocity field. Note: handling of oldTime faceU
    // fields happens in ddt scheme when calling ddtConsistentPhiCorr
    faceU = fvc::interpolate(U);

    // Interpolate original rAU on the faces
    const surfaceScalarField rAUf = fvc::interpolate(rAU);

    // Store previous iteration for the correct handling of under-relaxation
    phi.storePrevIter();

    // Calculate the ordinary part of the flux (H/A)
    phi = (faceU & mesh_.Sf());

    // STAGE 1: consistent ddt discretisation handling

    // Add consistent flux contribution due to ddt discretisation
    phi += fvc::ddtConsistentPhiCorr(faceU, U, rAUf);

    // Calculate aCoeff according to ddt discretisation
    aCoeff = 1.0 + ddtUEqn.A()*rAU;

    // STAGE 2: consistent under-relaxation handling

    // Get under-relaxation factor used in this iteration
    const dimensionedScalar alphaU = this->relaxFactor(U);

    // Helper variable
    const dimensionedScalar urfCoeff = (1.0 - alphaU)/alphaU;

    // Add consistent flux contribution due to possible under-relaxation
    phi += urfCoeff*fvc::interpolate(aCoeff)*phi.prevIter();

    // Add under-relaxation contribution to aCoeff
    aCoeff += urfCoeff*aCoeff;

    // Note: this is not 100% entirely consistent with the way FOAM handles
    // under-relaxation. In fvMatrix::relax() member function, the diagonal
    // dominance is first assured by limiting the diagonal and the matrix is
    // then relaxed. However, this is a very minor inconsistency, especially
    // since this term vanishes when non-linear iterations converge to a tight
    // tolerance (achieving steady state in SIMPLE or converging PISO/PIMPLE in
    // each time step). VV, 23/Dec/2016.

    // STAGE 3: scale the flux and correct it at the boundaries

    // Scale the flux
    phi /= fvc::interpolate(aCoeff);

    // Get necessary data at the boundaries
    const volVectorField::GeometricBoundaryField& Ub = U.boundaryField();
    const surfaceVectorField::GeometricBoundaryField& Sb =
        mesh_.Sf().boundaryField();

    surfaceScalarField::GeometricBoundaryField& phib = phi.boundaryField();

    // Correct flux at boundaries depending on patch type. Note: it is possible
    // that we missed some important boundary conditions that need special
    // considerations when calculating the flux. Maybe we can reorganise this in
    // future by relying on distinction between = and == operators for
    // fvsPatchFields. However, that would require serious changes at the
    // moment. VV, 23/Dec/2016.
    forAll (phib, patchI)
    {
        if (Ub[patchI].fixesValue())
        {
            // This is fixed value patch, flux needs to be recalculated
            // with respect to the boundary condition
            phib[patchI] == (Ub[patchI] & Sb[patchI]);
        }
        else if
        (
            isA<slipFvPatchVectorField>(Ub[patchI])
         || isA<symmetryFvPatchVectorField>(Ub[patchI])
        )
        {
            // This is slip or symmetry, flux needs to be zero
            phib[patchI] == 0.0;
        }
    }
}


void Foam::solutionControl::reconstructVelocity
(
    volVectorField& U,
    const fvVectorMatrix& ddtUEqn,
    const volScalarField& rAU,
    const volScalarField& p,
    const surfaceScalarField& phi
) const
{
    // Reconstruct the velocity using all the components from original equation
    U = 1.0/(1.0/rAU + ddtUEqn.A())*
        (
            U/rAU + ddtUEqn.H() - fvc::grad(p)
        );
    U.correctBoundaryConditions();

    // Update divergence free face velocity field, whose value will be used in
    // the next time step
    if (!faceUPtr_)
    {
        FatalErrorIn
        (
            "void solutionControl::reconstructVelocity"
            "\n("
            "\n    volVectorField& U,"
            "\n    const volVectorField& ddtUEqn,"
            "\n    const volScalarField& rAU,"
            "\n    const volScalarField& p"
            "\n) const"
        )   << "faceUPtr_ not calculated. Make sure you have called"
            << " calculateTimeConsistentFlux(...) before calling this function."
            << exit(FatalError);
    }
    surfaceVectorField& faceU = *faceUPtr_;

    // First interpolate the reconstructed velocity on the faces
    faceU = fvc::interpolate(U);

    // Replace the normal component with conservative flux

    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceVectorField rSf = Sf/magSqr(Sf);

    // Subtract interpolated normal component
    faceU -= (Sf & faceU)*rSf;

    // Now that the normal component is zero, add the normal component from
    // conservative flux
    faceU += phi*rSf;
}


const Foam::volScalarField& Foam::solutionControl::aCoeff() const
{
    if (!aCoeffPtr_)
    {
        FatalErrorIn
        (
            "const volScalarField& solutionControl::aCoeff() const"
        )   << "aCoeffPtr_ not calculated. Make sure you have called"
            << " calculateTimeConsistentFlux(...) before calling aCoeff()."
            << exit(FatalError);
    }

    return *aCoeffPtr_;
}


// ************************************************************************* //
