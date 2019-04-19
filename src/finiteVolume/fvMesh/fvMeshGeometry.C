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

#include "fvMesh.H"
#include "foamTime.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fvMesh::makeSf() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeSf() const : "
            << "assembling face areas"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (SfPtr_)
    {
        FatalErrorIn("fvMesh::makeSf() const")
            << "face areas already exist"
            << abort(FatalError);
    }

    SfPtr_ = new slicedSurfaceVectorField
    (
        IOobject
        (
            "S",
            pointsInstance(),
            meshSubDir,
            *this
        ),
        *this,
        dimArea,
        faceAreas()
    );

    // Boundary update.  Used in complex geometries, eg. immersed boundary
    // HJ, 29/Nov/2017
    slicedSurfaceVectorField& S = *SfPtr_;

    forAll (S.boundaryField(), patchI)
    {
        boundary()[patchI].makeSf(S);
    }
}


void fvMesh::makeMagSf() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeMagSf() const : "
            << "assembling mag face areas"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (magSfPtr_)
    {
        FatalErrorIn("void fvMesh::makeMagSf() const")
            << "mag face areas already exist"
            << abort(FatalError);
    }

    // Note: Added stabilisation for faces with exactly zero area.
    // These should be caught on mesh checking but at least this stops
    // the code from producing NaNs.  HJ, date deleted
    magSfPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "magSf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mag(Sf()) + dimensionedScalar("vs", dimArea, VSMALL)
    );

    // Note: boundary update not required, as magSf is calculated from Sf
    // HJ, 29/Nov/2017
}


void fvMesh::makeC() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeC() const : "
            << "assembling cell centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CPtr_)
    {
        FatalErrorIn("fvMesh::makeC() const")
            << "cell centres already exist"
            << abort(FatalError);
    }

    CPtr_ = new slicedVolVectorField
    (
        IOobject
        (
            "C",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimLength,
        cellCentres(),
        faceCentres(),
        false             // Do not preserve couples in geometry fields
    );

    // This piece of code is necessary for cyclic and cyclicGgi interfaces
    // using a separationOffset transform.
    // Those two interfaces will be using the method ::patchNeighbourField()
    // to evaluate the field C on the shadow patch. For cyclic and cyclicGgi
    // translational-only interfaces, the separationOffset transform is never
    // applied directly in ::patchNeighbourField() because this transform is
    // only pertinent for 3D coordinates, and the method ::patchNeighbourField()
    // does not discriminate the type of field it is operating on.
    // So, because the separationOffset transform is not applied, the evaluation
    // of a 3D position field like 'C' will always be wrong on the shadow
    // patches of translational cyclic and cyclicGgi interfaces.
    // For cyclic and cyclicGgi interfaces using a rotational transform, the
    // evaluation of the field C will be valid, but since we are only
    // interested in the patch face centers for these interfaces, we can
    // override those values as well.
    // See also:
    //    https://sourceforge.net/apps/mantisbt/openfoam-extend/view.php?id=42
    // MB, 12/Dec/2010
    //
    // Need to correct for cyclics transformation since absolute quantity.
    // Ok on processor patches since hold opposite cell centre (no
    // transformation)

    // Note: moved into virtual functions
    // HJ, 29/Nov/2017

    // Boundary update.  Used in complex geometries, eg. immersed boundary
    // HJ, 29/Nov/2017
    slicedVolVectorField& C = *CPtr_;

    forAll (C.boundaryField(), patchI)
    {
        boundary()[patchI].makeC(C);
    }
    // Note:
    // Functionality for cyclic and cyclicGgiFvPatch, which used to be here
    // with RTTI is moved into respective patched under virtual functions
    // HJ, 29/Nov/2017
}


void fvMesh::makeCf() const
{
    if (debug)
    {
        Info<< "void fvMesh::makeCf() const : "
            << "assembling face centres"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (CfPtr_)
    {
        FatalErrorIn("fvMesh::makeCf() const")
            << "face centres already exist"
            << abort(FatalError);
    }

    CfPtr_ = new slicedSurfaceVectorField
    (
        IOobject
        (
            "Cf",
            pointsInstance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        *this,
        dimLength,
        faceCentres(),
        false             // Do not preserve couples in geometry fields
    );

    // Boundary update.  Used in complex geometries, eg. immersed boundary
    // HJ, 29/Nov/2017
    slicedSurfaceVectorField& Cf = *CfPtr_;

    forAll (Cf.boundaryField(), patchI)
    {
        boundary()[patchI].makeCf(Cf);
    }
}


void fvMesh::makePhi() const
{
    if (debug)
    {
        InfoIn("void fvMesh::makePhi() const")
            << "Preparing mesh flux field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (phiPtr_)
    {
        FatalErrorIn("fvMesh::makePhi() const")
            << "flux field already exists"
            << abort(FatalError);
    }

    // Reading old time mesh motion flux if it exists and
    // creating zero current time mesh motion flux

    scalar t0 = time().value() - time().deltaT().value();

    IOobject meshPhiHeader
    (
        "meshPhi",
        time().timeName(t0),
        *this,
        IOobject::NO_READ
    );

    if (meshPhiHeader.headerOk())
    {
        if (debug)
        {
            InfoIn("void fvMesh::makePhi()")
                << "Reading mesh fluxes" << endl;
        }

        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                time().timeName(t0),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );

        phiPtr_->oldTime();

        (*phiPtr_) = dimensionedScalar("0", dimVolume/dimTime, 0.0);

        // This mesh is moving: set the motion to true
    }
    else
    {
        if (debug)
        {
            InfoIn("void fvMesh::makePhi()")
                << "Creating null mesh motion fluxes" << endl;
        }

        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimensionedScalar("0", dimVolume/dimTime, 0.0)
        );
    }
}


void fvMesh::updatePhi(const scalarField& sweptVols) const
{
    // Fill in mesh motion fluxes given swept volumes for all faces

    // Make sure V and V0 are constructed before the correction
    // HJ, 22/Dec/2017
    V();
    if (!V0Ptr_)
    {
        const_cast<fvMesh&>(*this).setV0();
    }

    if (!phiPtr_)
    {
        makePhi();
    }

    scalar rDeltaT = 1.0/time().deltaT().value();

    surfaceScalarField& phi = *phiPtr_;

    phi.internalField() = scalarField::subField(sweptVols, nInternalFaces());
    phi.internalField() *= rDeltaT;

    const fvPatchList& patches = boundary();

    // Calculate regular values first and then allow patches to update them
    // HJ, 15/Dec/2017
    forAll (patches, patchI)
    {
        phi.boundaryField()[patchI] = patches[patchI].patchSlice(sweptVols);
        phi.boundaryField()[patchI] *= rDeltaT;
    }

    // Boundary update.  Used in complex geometry updates, eg. immersed boundary
    // HJ, 29/Nov/2017
    forAll (phi.boundaryField(), patchI)
    {
        boundary()[patchI].updatePhi
        (
            *VPtr_,
            *V0Ptr_,
            phi
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volScalarField::DimensionedInternalField& fvMesh::V() const
{
    if (!VPtr_)
    {
        if (debug)
        {
            InfoIn
            (
                "const volScalarField::DimensionedInternalField& "
                "fvMesh::V() const"
            )   << "Calculating cell volumes." << endl;
        }

        VPtr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            *this,
            dimVolume,
            cellVolumes()
        );

        // Boundary update.  Used in complex geometries, eg. immersed boundary
        // HJ, 29/Nov/2017
        scalarField& V = *VPtr_;

        forAll (boundary(), patchI)
        {
            boundary()[patchI].makeV(V);
        }
    }

    return *VPtr_;
}


const volScalarField::DimensionedInternalField& fvMesh::V0() const
{
    if (!V0Ptr_)
    {
        FatalErrorIn("fvMesh::V0() const")
            << "V0 is not available"
            << abort(FatalError);
    }

    return *V0Ptr_;
}


DimensionedField<scalar, volMesh>& fvMesh::setV0()
{
    // Delete old volume and mesh motion fluxes.  setV0() must be followed by
    // another mesh motion.  HJ, 25/Feb/2009
    deleteDemandDrivenData(phiPtr_);
    deleteDemandDrivenData(V0Ptr_);

    if (debug)
    {
        InfoIn("DimensionedField<scalar, volMesh>& fvMesh::setV0()")
            << "Setting old cell volumes" << endl;
    }

    // Update time index
    curTimeIndex_ = time().timeIndex();

    V0Ptr_ = new DimensionedField<scalar, volMesh>
    (
        IOobject
        (
            "V0",
            time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        V()
    );

    return *V0Ptr_;
}


const volScalarField::DimensionedInternalField& fvMesh::V00() const
{
    if (!V00Ptr_)
    {
        V00Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V00",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            V0()
        );

        // If V00 is used then V0 should be stored for restart
        V0Ptr_->writeOpt() = IOobject::AUTO_WRITE;
    }

    return *V00Ptr_;
}


tmp<volScalarField::DimensionedInternalField> fvMesh::Vsc() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar tFrac =
        (
            ts.value() - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (tFrac < (1 - SMALL))
        {
            return V0() + tFrac*(V() - V0());
        }
        else
        {
            return V();
        }
    }
    else
    {
        return V();
    }
}


tmp<volScalarField::DimensionedInternalField> fvMesh::Vsc0() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar t0Frac =
        (
            (ts.value() - ts.deltaTValue())
          - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (t0Frac > SMALL)
        {
            return V0() + t0Frac*(V() - V0());
        }
        else
        {
            return V0();
        }
    }
    else
    {
        return V0();
    }
}


const surfaceVectorField& fvMesh::Sf() const
{
    if (!SfPtr_)
    {
        makeSf();
    }

    return *SfPtr_;
}


const surfaceScalarField& fvMesh::magSf() const
{
    if (!magSfPtr_)
    {
        makeMagSf();
    }

    return *magSfPtr_;
}


const volVectorField& fvMesh::C() const
{
    if (!CPtr_)
    {
        makeC();
    }

    return *CPtr_;
}


const surfaceVectorField& fvMesh::Cf() const
{
    if (!CfPtr_)
    {
        makeCf();
    }

    return *CfPtr_;
}


const surfaceScalarField& fvMesh::phi() const
{
    if (!phiPtr_)
    {
        // If making mesh motion fluxes from nothing, old volumes
        // must be recorded as well.  HJ, 10/Aug/2017
        const_cast<fvMesh&>(*this).setV0();
        makePhi();
    }

    // Set zero current time
    // mesh motion fluxes if the time has been incremented
    if (phiPtr_->timeIndex() != time().timeIndex())
    {
        phiPtr_->oldTime();

        if (debug)
        {
            InfoIn("const surfaceScalarField& fvMesh::phi() const")
                << "Resetting mesh motion fluxes to zero" << endl;
        }

        (*phiPtr_) = dimensionedScalar("0", dimVolume/dimTime, 0.0);
    }

    return *phiPtr_;
}


surfaceScalarField& fvMesh::setPhi()
{
    if (!phiPtr_)
    {
        makePhi();
    }

    return *phiPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
