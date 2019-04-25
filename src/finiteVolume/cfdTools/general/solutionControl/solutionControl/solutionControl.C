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

#include "solutionControl.H"
#include "fieldTypes.H"
#include "VectorNFieldTypes.H"
#include "fvc.H"
#include "inletOutletFvPatchFields.H"
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

    storePrevIter<vector2>();
    storePrevIter<vector4>();
    storePrevIter<vector6>();
    storePrevIter<vector8>();

    storePrevIter<sphericalTensor2>();
    storePrevIter<sphericalTensor4>();
    storePrevIter<sphericalTensor6>();
    storePrevIter<sphericalTensor8>();

    storePrevIter<tensor2>();
    storePrevIter<tensor4>();
    storePrevIter<tensor6>();
    storePrevIter<tensor8>();
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
        const List<BlockSolverPerformance<Type> > sp(data);
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

    maxTypeResidual<vector2>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<vector4>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<vector6>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<vector8>(fieldName, data, firstRes, lastRes);

    maxTypeResidual<sphericalTensor2>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<sphericalTensor4>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<sphericalTensor6>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<sphericalTensor8>(fieldName, data, firstRes, lastRes);

    maxTypeResidual<tensor2>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<tensor4>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<tensor6>(fieldName, data, firstRes, lastRes);
    maxTypeResidual<tensor8>(fieldName, data, firstRes, lastRes);

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


void Foam::solutionControl::addDdtFluxContribution
(
    surfaceScalarField& phi,
    surfaceScalarField& aCoeff,
    const surfaceVectorField& faceU,
    const volVectorField& U,
    const surfaceScalarField& rAUf,
    const fvVectorMatrix& ddtUEqn
) const
{
    // Add ddt scheme dependent flux contribution. Here: phi = H/A & Sf.
    phi += fvc::ddtConsistentPhiCorr(faceU, U, rAUf);

    // Add ddt scheme dependent contribution to diffusion scale coeff. Here:
    // aCoeff = 1.
    aCoeff += fvc::interpolate(ddtUEqn.A())*rAUf;
}


void Foam::solutionControl::addUnderRelaxationFluxContribution
(
    surfaceScalarField& phi,
    surfaceScalarField& aCoeff,
    const volVectorField& U
) const
{
    // Get under-relaxation factor used in this iteration
    const dimensionedScalar alphaU(this->relaxFactor(U));

    // Return if the under-relaxation factor is >= 1 and =< 0
    if ((alphaU.value() > 1.0 - SMALL) || (alphaU.value() < SMALL))
    {
        return;
    }

    const dimensionedScalar urfCoeff((1.0 - alphaU)/alphaU);

    // Add under-relaxation dependent flux contribution
    phi += urfCoeff*aCoeff*phi.prevIter();

    // Add under-relaxation dependent contribution to diffusion scale coeff
    aCoeff += urfCoeff*aCoeff;
}


void Foam::solutionControl::correctBoundaryFlux
(
    surfaceScalarField& phi,
    const volVectorField& U
) const
{
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
    // moment (e.g. using phi = rUA*UEqn.H() as in most solvers would not update
    // the boundary values correctly for fixed fvsPatchFields). VV, 20/Apr/2017.
    forAll (phib, patchI)
    {
        // Get patch field
        const fvPatchVectorField& Up = Ub[patchI];

        if
        (
            Up.fixesValue()
         && !isA<inletOutletFvPatchVectorField>(Up)
        )
        {
            // This is fixed value patch, flux needs to be recalculated
            // with respect to the boundary condition
            phib[patchI] == (Up & Sb[patchI]);
        }
        else if
        (
            isA<slipFvPatchVectorField>(Up)
         || isA<symmetryFvPatchVectorField>(Up)
        )
        {
            // This is slip or symmetry, flux needs to be zero
            phib[patchI] == 0.0;
        }
    }
}


void Foam::solutionControl::correctBoundaryFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const surfaceScalarField& meshPhi
) const
{
    // Get necessary data at the boundaries
    const volVectorField::GeometricBoundaryField& Ub = U.boundaryField();
    const surfaceVectorField::GeometricBoundaryField& Sb =
        mesh_.Sf().boundaryField();

    surfaceScalarField::GeometricBoundaryField& phib = phi.boundaryField();

    const surfaceScalarField::GeometricBoundaryField& meshPhib =
        meshPhi.boundaryField();

    // Correct flux at boundaries depending on patch type. Note: it is possible
    // that we missed some important boundary conditions that need special
    // considerations when calculating the flux. Maybe we can reorganise this in
    // future by relying on distinction between = and == operators for
    // fvsPatchFields. However, that would require serious changes at the
    // moment (e.g. using phi = rUA*UEqn.H() as in most solvers would not update
    // the boundary values correctly for fixed fvsPatchFields). VV, 20/Apr/2017.
    forAll (phib, patchI)
    {
        // Get patch field
        const fvPatchVectorField& Up = Ub[patchI];

        if
        (
            Up.fixesValue()
         && !isA<inletOutletFvPatchVectorField>(Up)
        )
        {
            // This is fixed value patch, flux needs to be recalculated
            // with respect to the boundary condition
            phib[patchI] == (Up & Sb[patchI]) - meshPhib[patchI];
        }
        else if
        (
            isA<slipFvPatchVectorField>(Up)
         || isA<symmetryFvPatchVectorField>(Up)
        )
        {
            // This is slip or symmetry, flux needs to be zero
            phib[patchI] == 0.0;
        }
    }
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
    aCoeffPtrs_(),
    faceUPtrs_(),
    indices_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionControl::~solutionControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solutionControl::calcTransientConsistentFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& rAU,
    const fvVectorMatrix& ddtUEqn
) const
{
    // Store necessary data for this velocity field
    const word& UName = U.name();

    // Check whether the fields are present in the list
    if (!indices_.found(UName))
    {
        // Get current index as size of the indices list (before insertion)
        const label i = indices_.size();

        // Insert the index into the hash table
        indices_.insert(UName, i);

        // Extend lists
        aCoeffPtrs_.resize(indices_.size());
        faceUPtrs_.resize(indices_.size());

        // Double check whether the fields have been already set
        if (!aCoeffPtrs_.set(i) && !faceUPtrs_.set(i))
        {
            aCoeffPtrs_.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "aCoeff." + UName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );

            faceUPtrs_.set
            (
                i,
                new surfaceVectorField
                (
                    IOobject
                    (
                        "faceU." + UName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("zero", dimVelocity, vector::zero)
                )
            );
        }
        else if (!aCoeffPtrs_.set(i) || !faceUPtrs_.set(i))
        {
            FatalErrorIn
            (
                "void solutionControl::calcTransientConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const volScalarField& rAU,"
                "\n    const fvVectorMatrix& ddtUEqn"
                "\n)"
            )   << "Either aCoeff or faceU is allocated for field " << UName
                << " while the other is not." << nl
                << " This must not happen in transient simulation. Make sure"
                << " that functions aiding consistency are called in the right"
                << " order (first flux and then velocity reconstruction)."
                << exit(FatalError);
        }
    }
    else
    {
        // Index has been set for this field, so the fields must be there as
        // well. Check and report an error if they are not allocated
        if (!aCoeffPtrs_.set(indices_[UName]))
        {
            FatalErrorIn
            (
                "void solutionControl::calcTransientConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const volScalarField& rAU,"
                "\n    const fvVectorMatrix& ddtUEqn"
                "\n)"
            )   << "Index is set, but the aCoeff field is not allocated for "
                << UName << "." << nl
                << "This should not happen for transient simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
        else if (!faceUPtrs_.set(indices_[UName]))
        {
            FatalErrorIn
            (
                "void solutionControl::calcTransientConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const volScalarField& rAU,"
                "\n    const fvVectorMatrix& ddtUEqn"
                "\n)"
            )   << "Index is set, but the faceU field is not allocated for "
                << UName << "." << nl
                << "This should not happen for transient simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
    }

    // Algorithm:
    // 1. Update flux and aCoeff due to ddt discretisation
    // 2. Update flux and aCoeff due to under-relaxation
    // 3. Scale the flux with aCoeff, making sure that flux at fixed boundaries
    //    remains consistent

    // Get index from the hash table
    const label i = indices_[UName];

    // Get fields that will be updated
    surfaceScalarField& aCoeff = aCoeffPtrs_[i];
    surfaceVectorField& faceU = faceUPtrs_[i];

    // Update face interpolated velocity field. Note: handling of oldTime faceU
    // fields happens in ddt scheme when calling ddtConsistentPhiCorr inside
    // addDdtFluxContribution
    faceU = fvc::interpolate(U);

    // Interpolate original rAU on the faces
    const surfaceScalarField rAUf = fvc::interpolate(rAU);

    // Store previous iteration for the correct handling of under-relaxation
    phi.storePrevIter();

    // Calculate the ordinary part of the flux (H/A)
    phi = (faceU & mesh_.Sf());

    // Initialize aCoeff to 1
    aCoeff = dimensionedScalar("one", dimless, 1.0);

    // STAGE 1: consistent ddt discretisation handling
    addDdtFluxContribution(phi, aCoeff, faceU, U, rAUf, ddtUEqn);

    // STAGE 2: consistent under-relaxation handling
    addUnderRelaxationFluxContribution(phi, aCoeff, U);

    // STAGE 3: scale the flux and correct it at the boundaries
    phi /= aCoeff;
    correctBoundaryFlux(phi, U);
}


void Foam::solutionControl::calcSteadyConsistentFlux
(
    surfaceScalarField& phi,
    const volVectorField& U
) const
{
    // Store necessary data for this velocity field
    const word& UName = U.name();

    // Check whether the fields are present in the list
    if (!indices_.found(UName))
    {
        // Get current index as size of the indices list (before insertion)
        const label i = indices_.size();

        // Insert the index into the hash table
        indices_.insert(UName, i);

        // Extend list (only aCoeff list because faceU does not need to be
        // stored for steady state simulations)
        aCoeffPtrs_.resize(indices_.size());

        // Double check whether the aCoeff field has been already set. Note:
        // faceU does not need to be stored for steady state
        if (!aCoeffPtrs_.set(i) && faceUPtrs_.empty())
        {
            aCoeffPtrs_.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "aCoeff." + UName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );
        }
        else if (!aCoeffPtrs_.set(i) || !faceUPtrs_.empty())
        {
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U"
                "\n)"
            )   << "Either aCoeff or faceU is allocated for field " << UName
                << " while the other is not."
                << "This must not happen in transient simulation. Make sure"
                << " that functions aiding consistency are called in the right"
                << " order (first flux and then velocity reconstruction)."
                << exit(FatalError);
        }
    }
    else
    {
        // Index has been set for this field. Check the state of fields
        if (!aCoeffPtrs_.set(indices_[UName]))
        {
            // aCoeff field should be allocated
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U"
                "\n)"
            )   << "Index is set, but the aCoeff field is not allocated for "
                << UName << "." << nl
                << "This should not happen for steady state simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
        else if (!faceUPtrs_.empty())
        {
            // faceU field shouldn't be allocated
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U"
                "\n)"
            )   << "Index is set, but the faceU field is allocated for "
                << UName << "." << nl
                << "This should not happen for steady state simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
    }

    // Algorithm:
    // 1. Update flux and aCoeff due to under-relaxation
    // 2. Scale the flux with aCoeff, making sure that flux at fixed boundaries
    //    remains consistent

    // Get index from the hash table
    const label i = indices_[UName];

    // Get aCoeff field. Note: no need to check whether the entry has been found
    // since we have just inserted it above
    surfaceScalarField& aCoeff = aCoeffPtrs_[i];

    // Store previous iteration for the correct handling of under-relaxation
    phi.storePrevIter();

    // Calculate the ordinary part of the flux (H/A)
    phi = (fvc::interpolate(U) & mesh_.Sf());

    // Initialize aCoeff to 1
    aCoeff = dimensionedScalar("one", dimless, 1.0);

    // STAGE 1: consistent under-relaxation handling
    addUnderRelaxationFluxContribution(phi, aCoeff, U);

    // STAGE 2: scale the flux and correct it at the boundaries
    phi /= aCoeff;
    correctBoundaryFlux(phi, U);
}


void Foam::solutionControl::calcSteadyConsistentFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const surfaceScalarField& meshPhi
) const
{
    // Store necessary data for this velocity field
    const word& UName = U.name();

    // Check whether the fields are present in the list
    if (!indices_.found(UName))
    {
        // Get current index as size of the indices list (before insertion)
        const label i = indices_.size();

        // Insert the index into the hash table
        indices_.insert(UName, i);

        // Extend list (only aCoeff list because faceU does not need to be
        // stored for steady state simulations)
        aCoeffPtrs_.resize(indices_.size());

        // Double check whether the aCoeff field has been already set. Note:
        // faceU does not need to be stored for steady state
        if (!aCoeffPtrs_.set(i) && faceUPtrs_.empty())
        {
            aCoeffPtrs_.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "aCoeff." + UName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );
        }
        else if (!aCoeffPtrs_.set(i) || !faceUPtrs_.empty())
        {
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U"
                "\n)"
            )   << "Either aCoeff or faceU is allocated for field " << UName
                << " while the other is not."
                << "This must not happen in transient simulation. Make sure"
                << " that functions aiding consistency are called in the right"
                << " order (first flux and then velocity reconstruction)."
                << exit(FatalError);
        }
    }
    else
    {
        // Index has been set for this field. Check the state of fields
        if (!aCoeffPtrs_.set(indices_[UName]))
        {
            // aCoeff field should be allocated
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U"
                "\n)"
            )   << "Index is set, but the aCoeff field is not allocated for "
                << UName << "." << nl
                << "This should not happen for steady state simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
        else if (!faceUPtrs_.empty())
        {
            // faceU field shouldn't be allocated
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U"
                "\n)"
            )   << "Index is set, but the faceU field is allocated for "
                << UName << "." << nl
                << "This should not happen for steady state simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
    }

    // Algorithm:
    // 1. Update flux and aCoeff due to under-relaxation
    // 2. Scale the flux with aCoeff, making sure that flux at fixed boundaries
    //    remains consistent

    // Get index from the hash table
    const label i = indices_[UName];

    // Get aCoeff field. Note: no need to check whether the entry has been found
    // since we have just inserted it above
    surfaceScalarField& aCoeff = aCoeffPtrs_[i];

    // Store previous iteration for the correct handling of under-relaxation
    phi.storePrevIter();

    // Calculate the ordinary part of the flux (H/A)
    phi = (fvc::interpolate(U) & mesh_.Sf()) - meshPhi;

//    surfaceScalarField meshPhi("meshPhi", phi);
//    mrfZones.relativeFlux(phi);

//    meshPhi -= phi;

    // Initialize aCoeff to 1
    aCoeff = dimensionedScalar("one", dimless, 1.0);

    // STAGE 1: consistent under-relaxation handling
    addUnderRelaxationFluxContribution(phi, aCoeff, U);

    // STAGE 2: scale the flux and correct it at the boundaries
    phi /= aCoeff;
    correctBoundaryFlux(phi, U, meshPhi);
}


void Foam::solutionControl::calcSteadyMRFConsistentFlux
(
    surfaceScalarField& phi,
    const volVectorField& U,
    const MRFZones& mrfZones
) const
{
    // Store necessary data for this velocity field
    const word& UName = U.name();

    // Check whether the fields are present in the list
    if (!indices_.found(UName))
    {
        // Get current index as size of the indices list (before insertion)
        const label i = indices_.size();

        // Insert the index into the hash table
        indices_.insert(UName, i);

        // Extend list (only aCoeff list because faceU does not need to be
        // stored for steady state simulations)
        aCoeffPtrs_.resize(indices_.size());

        // Double check whether the aCoeff field has been already set. Note:
        // faceU does not need to be stored for steady state
        if (!aCoeffPtrs_.set(i) && faceUPtrs_.empty())
        {
            aCoeffPtrs_.set
            (
                i,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "aCoeff." + UName,
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );
        }
        else if (!aCoeffPtrs_.set(i) || !faceUPtrs_.empty())
        {
            FatalErrorIn
            (
                "void solutionControl::calcSteadyMRFConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const MRFZones& mrfZones"
                "\n)"
            )   << "Either aCoeff or faceU is allocated for field " << UName
                << " while the other is not."
                << "This must not happen in transient simulation. Make sure"
                << " that functions aiding consistency are called in the right"
                << " order (first flux and then velocity reconstruction)."
                << exit(FatalError);
        }
    }
    else
    {
        // Index has been set for this field. Check the state of fields
        if (!aCoeffPtrs_.set(indices_[UName]))
        {
            // aCoeff field should be allocated
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const MRFZones& mrfZones"
                "\n)"
            )   << "Index is set, but the aCoeff field is not allocated for "
                << UName << "." << nl
                << "This should not happen for steady state simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
        else if (!faceUPtrs_.empty())
        {
            // faceU field shouldn't be allocated
            FatalErrorIn
            (
                "void solutionControl::calcSteadyConsistentFlux"
                "\n("
                "\n    surfaceScalarField& phi,"
                "\n    const volVectorField& U,"
                "\n    const MRFZones& mrfZones"
                "\n)"
            )   << "Index is set, but the faceU field is allocated for "
                << UName << "." << nl
                << "This should not happen for steady state simulation." << nl
                << "Something went wrong."
                << exit(FatalError);
        }
    }

    // Algorithm:
    // 1. Update flux and aCoeff due to under-relaxation
    // 2. Scale the flux with aCoeff, making sure that flux at fixed boundaries
    //    remains consistent

    // Get index from the hash table
    const label i = indices_[UName];

    // Get aCoeff field. Note: no need to check whether the entry has been found
    // since we have just inserted it above
    surfaceScalarField& aCoeff = aCoeffPtrs_[i];

    // Store previous iteration for the correct handling of under-relaxation
    phi.storePrevIter();

    // Calculate the ordinary part of the flux (H/A)
    phi = (fvc::interpolate(U) & mesh_.Sf());

    surfaceScalarField meshPhi("meshPhi", phi);
    mrfZones.relativeFlux(phi);

    meshPhi -= phi;

    // Initialize aCoeff to 1
    aCoeff = dimensionedScalar("one", dimless, 1.0);

    // STAGE 1: consistent under-relaxation handling
    addUnderRelaxationFluxContribution(phi, aCoeff, U);

    // STAGE 2: scale the flux and correct it at the boundaries
    phi /= aCoeff;
    correctBoundaryFlux(phi, U, meshPhi);
}


void Foam::solutionControl::reconstructTransientVelocity
(
    volVectorField& U,
    surfaceScalarField& phi,
    const fvVectorMatrix& ddtUEqn,
    const volScalarField& rAU,
    const volScalarField& p
) const
{
    // Reconstruct the velocity using all the components from original equation
    U = 1.0/(1.0/rAU + ddtUEqn.A())*
        (
            U/rAU + ddtUEqn.H() - fvc::grad(p)
        );

    // Get name and the corresponding index
    const word& UName = U.name();
    const label i = indices_[UName];

    // Update divergence free face velocity field, whose value will be used in
    // the next time step. Note that we need to have absolute flux here.
    if (!faceUPtrs_.set(i))
    {
        FatalErrorIn
        (
            "void solutionControl::reconstructTransientVelocity"
            "\n("
            "\n    volVectorField& U,"
            "\n    const surfaceScalarField& phi,"
            "\n    const volVectorField& ddtUEqn,"
            "\n    const volScalarField& rAU,"
            "\n    const volScalarField& p"
            "\n) const"
        )   << "faceU not calculated for field " << UName
            << ". Make sure you have called"
            << " calcTransientConsistentFlux(...) before calling this function."
            << exit(FatalError);
    }

    // Get faceU field. Note: no need to check whether the entry has been found
    // since we have just checked
    surfaceVectorField& faceU = faceUPtrs_[i];

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

    // If the mesh is moving, flux needs to be relative before boundary
    // conditions for velocity are corrected. VV and IG, 4/Jan/2016.
    fvc::makeRelative(phi, U);

    // Correct boundary conditions with relative flux
    U.correctBoundaryConditions();
}


void Foam::solutionControl::reconstructSteadyVelocity
(
    volVectorField& U,
    const volScalarField& rAU,
    const volScalarField& p
) const
{
    // Reconstruct the velocity field
    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();

    // Note: no need to store and update faceU field for steady state run
}


const Foam::surfaceScalarField& Foam::solutionControl::aCoeff
(
    const word& UName
) const
{
    // Get corresponding index
    const label i = indices_[UName];

    if (!aCoeffPtrs_.set(i))
    {
        FatalErrorIn
        (
            "const surfaceScalarField& solutionControl::aCoeff"
            "(const word& UName) const"
        )   << "aCoeff not calculated for field " << UName
            << ". Make sure you have called"
            << " calcTransientConsistentFlux(...) or "
            << " calcSteadyConsistentFlux(...) before calling aCoeff()."
            << exit(FatalError);
    }

    return aCoeffPtrs_[i];
}


// ************************************************************************* //
