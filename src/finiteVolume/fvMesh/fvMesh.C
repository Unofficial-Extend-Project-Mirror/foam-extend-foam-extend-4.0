/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "fvMeshLduAddressing.H"
#include "emptyPolyPatch.H"
#include "mapPolyMesh.H"
#include "MapFvFields.H"
#include "fvMeshMapper.H"
#include "mapClouds.H"
#include "meshObjectBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fvMesh, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void fvMesh::clearGeomNotOldVol()
{
    deleteDemandDrivenData(VPtr_);

    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfPtr_);
}


void fvMesh::clearGeom()
{
    clearGeomNotOldVol();

    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);

    // Mesh motion flux cannot be deleted here because the old-time flux
    // needs to be saved.
}


void fvMesh::clearAddressing()
{
    deleteDemandDrivenData(lduPtr_);
}


void fvMesh::clearOut()
{
    clearGeom();
    surfaceInterpolation::clearOut();

    clearAddressing();

    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);

    polyMesh::clearOut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fvMesh::fvMesh(const IOobject& io)
:
    polyMesh(io),
    surfaceInterpolation(*this),
    boundary_(*this, boundaryMesh()),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    VPtr_(NULL),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from IOobject"
            << endl;
    }

    // Check the existance of the cell volumes and read if present
    // and set the storage of V00
    if (file(time().timePath()/"V0"))
    {
        if (debug)
        {
            Info<< "Reading old cell volumes" << endl;
        }

        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V0",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            *this
        );

        V00();
    }

    // Check the existance of the mesh fluxes, read if present and set the 
    // mesh to be moving
    if (file(time().timePath()/"meshPhi"))
    {
        if (debug)
        {
            Info<< "Reading motion fluxes" << endl;
        }

        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );

        // The mesh is now considered moving so the old-time cell volumes
        // will be required for the time derivatives so if they haven't been
        // read initialise to the current cell volumes
        if (!V0Ptr_)
        {
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "V0",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                V()
            );
        }

        moving(true);
    }
}


fvMesh::fvMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const labelList& allOwner,
    const labelList& allNeighbour,
    const bool syncPar
)
:
    polyMesh(io, points, faces, allOwner, allNeighbour, syncPar),
    surfaceInterpolation(*this),
    boundary_(*this),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    VPtr_(NULL),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from components"
            << endl;
    }
}


fvMesh::fvMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const bool syncPar
)
:
    polyMesh(io, points, faces, cells, syncPar),
    surfaceInterpolation(*this),
    boundary_(*this),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    VPtr_(NULL),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from components"
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Helper function for construction from pieces
void fvMesh::addFvPatches(const List<polyPatch*> & p, const bool validBoundary)
{
    if (boundary().size() > 0)
    {
        FatalErrorIn
        (
            "fvMesh::addFvPatches(const List<polyPatch*>&, const bool)"
        )   << " boundary already exists"
            << abort(FatalError);
    }

    // first add polyPatches
    addPatches(p, validBoundary);
    boundary_.addPatches(boundaryMesh());
}


void fvMesh::removeFvBoundary()
{
    if (debug)
    {
        Info<< "void fvMesh::removeFvBoundary(): "
            << "Removing boundary patches."
            << endl;
    }

    // Remove fvBoundaryMesh data first.
    boundary_.clear();
    boundary_.setSize(0);
    polyMesh::removeBoundary();

    clearOut();
}


polyMesh::readUpdateState fvMesh::readUpdate()
{
    if (debug)
    {
        Info<< "polyMesh::readUpdateState fvMesh::readUpdate() : "
            << "Updating fvMesh.  ";
    }

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        if (debug)
        {
            Info << "Boundary and topological update" << endl;
        }

        boundary_.readUpdate(boundaryMesh());

        clearOut();
    }
    else if (state == polyMesh::TOPO_CHANGE)
    {
        if (debug)
        {
            Info << "Topological update" << endl;
        }

        clearOut();
    }
    else if (state == polyMesh::POINTS_MOVED)
    {
        if (debug)
        {
            Info << "Point motion update" << endl;
        }

        clearGeom();
    }
    else
    {
        if (debug)
        {
            Info << "No update" << endl;
        }
    }

    return state;
}


const fvBoundaryMesh& fvMesh::boundary() const
{
    return boundary_;
}


const lduAddressing& fvMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new fvMeshLduAddressing(*this);
    }

    return *lduPtr_;
}


void fvMesh::mapFields(const mapPolyMesh& meshMap)
{
    if (debug)
    {
        Info<< "void fvMesh::mapFields(const mapPolyMesh& meshMap): "
            << "Mapping fv fields."
            << endl;
    }

    // Create a mapper
    const fvMeshMapper mapper(*this, meshMap);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar, fvPatchField, fvMeshMapper, volMesh>(mapper);
    MapGeometricFields<vector, fvPatchField, fvMeshMapper, volMesh>(mapper);
    MapGeometricFields<sphericalTensor, fvPatchField, fvMeshMapper, volMesh>
        (mapper);
    MapGeometricFields<symmTensor, fvPatchField, fvMeshMapper, volMesh>(mapper);
    MapGeometricFields<tensor, fvPatchField, fvMeshMapper, volMesh>(mapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar, fvsPatchField, fvMeshMapper, surfaceMesh>
        (mapper);
    MapGeometricFields<vector, fvsPatchField, fvMeshMapper, surfaceMesh>
        (mapper);
    MapGeometricFields
        <sphericalTensor, fvsPatchField, fvMeshMapper, surfaceMesh>(mapper);
    MapGeometricFields<symmTensor, fvsPatchField, fvMeshMapper, surfaceMesh>
        (mapper);
    MapGeometricFields<tensor, fvsPatchField, fvMeshMapper, surfaceMesh>
        (mapper);

    // Map all the clouds in the objectRegistry
    mapClouds(*this, meshMap);
}


void fvMesh::mapOldVolumes(const mapPolyMesh& meshMap)
{
    const labelList& cellMap = meshMap.cellMap();

    // Map the old volume. Just map to new cell labels.
    if (V0Ptr_)
    {
        if (debug)
        {
            InfoIn("void fvMesh::mapOldVolumes(const mapPolyMesh& meshMap)")
                << "Mapping old cell volumes." << endl;
        }

        scalarField& V0 = *V0Ptr_;

        scalarField savedV0(V0);
        V0.setSize(nCells());

        forAll(V0, i)
        {
            if (cellMap[i] > -1)
            {
                V0[i] = savedV0[cellMap[i]];
            }
            else
            {
                V0[i] = 0.0;
            }
        }
    }

    // Map the old-old volume. Just map to new cell labels.
    if (V00Ptr_)
    {
        if (debug)
        {
            InfoIn("void fvMesh::mapOldVolumes(const mapPolyMesh& meshMap)")
                << "Mapping old-old cell volumes." << endl;
        }

        scalarField& V00 = *V00Ptr_;

        scalarField savedV00(V00);
        V00.setSize(nCells());

        forAll(V00, i)
        {
            if (cellMap[i] > -1)
            {
                V00[i] = savedV00[cellMap[i]];
            }
            else
            {
                V00[i] = 0.0;
            }
        }
    }
}


void fvMesh::updateMesh(const mapPolyMesh& mpm)
{
    // Update polyMesh. This needs to keep volume existent!
    polyMesh::updateMesh(mpm);

    surfaceInterpolation::clearOut();
    clearGeomNotOldVol();

    // Map all fields
    mapFields(mpm);

    // Map old-volumes
    mapOldVolumes(mpm);

    clearAddressing();

    // handleMorph() should also clear out the surfaceInterpolation.
    // This is a temporary solution
    surfaceInterpolation::movePoints();

    meshObjectBase::allUpdateTopology<fvMesh>(*this, mpm);
}


void fvMesh::syncUpdateMesh()
{
    // Update polyMesh. This needs to keep volume existent!
    polyMesh::syncUpdateMesh();

    // Not sure how much clean-up is needed here.  HJ, 27/Nov/2009

    surfaceInterpolation::clearOut();
    clearGeomNotOldVol();

    clearAddressing();

    // handleMorph() should also clear out the surfaceInterpolation.
    // This is a temporary solution
    surfaceInterpolation::movePoints();

    // Instantiate a dummy mapPolyMesh
    autoPtr<mapPolyMesh> mapPtr(new mapPolyMesh(*this));

    meshObjectBase::allUpdateTopology<fvMesh>(*this, mapPtr());
}


tmp<scalarField> fvMesh::movePoints(const pointField& p)
{
    // Grab old time volumes if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        if (V00Ptr_ && V0Ptr_)
        {
            if (debug)
            {
                InfoIn("void fvMesh::movePoints(const mapPolyMesh& meshMap)")
                    << "Grabbing old-old cell volumes." << endl;
            }

            *V00Ptr_ = *V0Ptr_;
        }

        if (V0Ptr_)
        {
            if (debug)
            {
                InfoIn("void fvMesh::movePoints(const mapPolyMesh& meshMap)")
                    << "Grabbing old cell volumes." << endl;
            }

            *V0Ptr_ = V();
        }
        else
        {
            if (debug)
            {
                InfoIn("void fvMesh::movePoints(const mapPolyMesh& meshMap)")
                    << "Creating old cell volumes." << endl;
            }

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
        }

        curTimeIndex_ = time().timeIndex();
    }


    // Delete out of date geometrical information
    clearGeomNotOldVol();


    if (!phiPtr_)
    {
        // Create mesh motion flux
        if (debug)
        {
            InfoIn("tmp<scalarField> fvMesh::movePoints(const pointField& p)")
                << "Creating new mesh motion fluxes" << endl;
        }

        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimVolume/dimTime
        );
    }
    else
    {
        // Grab old time mesh motion fluxes if the time has been incremented
        if (phiPtr_->timeIndex() < time().timeIndex())
        {
            phiPtr_->oldTime();
        }
    }

    // Move the polyMesh and set the mesh motion fluxes to the swept-volumes
    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);

    updatePhi(tsweptVols());

    boundary_.movePoints();
    surfaceInterpolation::movePoints();
    meshObjectBase::allMovePoints<fvMesh>(*this);

    return tsweptVols;
}


bool fvMesh::writeObjects
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return polyMesh::writeObject(fmt, ver, cmp);
}


// Write mesh using IO settings from the time
bool fvMesh::write() const
{
    return polyMesh::write();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool fvMesh::operator!=(const fvMesh& bm) const
{
    return &bm != this;
}


bool fvMesh::operator==(const fvMesh& bm) const
{
    return &bm == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
