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

#include "polyMesh.H"
#include "Time.H"
#include "cellIOList.H"
#include "wedgePolyPatch.H"
#include "emptyPolyPatch.H"
#include "globalMeshData.H"
#include "processorPolyPatch.H"
#include "OSspecific.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMesh, 0);
}

Foam::word Foam::polyMesh::defaultRegion = "region0";
Foam::word Foam::polyMesh::meshSubDir = "polyMesh";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyMesh::calcDirections() const
{
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        directions_[cmpt] = 1;
    }

    label nEmptyPatches = 0;

    vector dirVec = vector::zero;

    forAll(boundaryMesh(), patchi)
    {
        if (isA<emptyPolyPatch>(boundaryMesh()[patchi]))
        {
            if (boundaryMesh()[patchi].size())
            {
                nEmptyPatches++;
                dirVec += sum(cmptMag(boundaryMesh()[patchi].faceAreas()));
            }
        }
    }

    if (nEmptyPatches)
    {
        reduce(dirVec, sumOp<vector>());

        dirVec /= mag(dirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (dirVec[cmpt] > 1e-6)
            {
                directions_[cmpt] = -1;
            }
            else
            {
                directions_[cmpt] = 1;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMesh::polyMesh(const IOobject& io)
:
    objectRegistry(io),
    primitiveMesh(),
    allPoints_
    (
        IOobject
        (
            "points",
            time().findInstance(meshDir(), "points"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    // To be re-sliced later.  HJ, 19/oct/2008
    points_(allPoints_, allPoints_.size()),
    allFaces_
    (
        IOobject
        (
            "faces",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    // To be re-sliced later.  HJ, 19/oct/2008
    faces_(allFaces_, allFaces_.size()),
    owner_
    (
        IOobject
        (
            "owner",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            time().findInstance(meshDir(), "faces"),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            time().findInstance(meshDir(), "boundary"),
            meshSubDir,
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    bounds_(allPoints_),
    directions_(Vector<label>::zero),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            time().findInstance
            (
                meshDir(),
                "pointZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            time().findInstance
            (
                meshDir(),
                "faceZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            time().findInstance
            (
                meshDir(),
                "cellZones",
                IOobject::READ_IF_PRESENT
            ),
            meshSubDir,
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    changing_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldAllPointsPtr_(NULL),
    oldPointsPtr_(NULL)
{
    if (exists(owner_.objectPath()))
    {
        initMesh();
    }
    else
    {
        cellIOList c
        (
            IOobject
            (
                "cells",
                // Find the cells file on the basis of the faces file
                // HJ, 8/Jul/2009
//                 time().findInstance(meshDir(), "cells"),
                time().findInstance(meshDir(), "faces"),
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );


        // Set the primitive mesh
        initMesh(c);

        owner_.write();
        neighbour_.write();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    // Warn if global empty mesh (constructs globalData!)
    if (globalData().nTotalPoints() == 0)
    {
        WarningIn("polyMesh(const IOobject&)")
            << "no points in mesh" << endl;
    }
    if (globalData().nTotalCells() == 0)
    {
        WarningIn("polyMesh(const IOobject&)")
            << "no cells in mesh" << endl;
    }
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const labelList& owner,
    const labelList& neighbour,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    allPoints_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    // To be re-sliced later.  HJ, 19/oct/2008
    points_(allPoints_, allPoints_.size()),
    allFaces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        faces
    ),
    // To be re-sliced later.  HJ, 19/oct/2008
    faces_(allFaces_, allFaces_.size()),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        owner
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        neighbour
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        0
    ),
    bounds_(allPoints_, syncPar),
    directions_(Vector<label>::zero),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    changing_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldAllPointsPtr_(NULL),
    oldPointsPtr_(NULL)
{
    // Check if the faces and cells are valid
    forAll (allFaces_, faceI)
    {
        const face& curFace = allFaces_[faceI];

        if (min(curFace) < 0 || max(curFace) > allPoints_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject& io,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const cellList& cells\n"
                ")\n"
            )   << "Face " << faceI << "contains vertex labels out of range: "
                << curFace << " Max point index = " << allPoints_.size()
                << abort(FatalError);
        }
    }

    // Set the primitive mesh
    initMesh();
}


Foam::polyMesh::polyMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const cellList& cells,
    const bool syncPar
)
:
    objectRegistry(io),
    primitiveMesh(),
    allPoints_
    (
        IOobject
        (
            "points",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        points
    ),
    // To be re-sliced later.  HJ, 19/oct/2008
    points_(allPoints_, allPoints_.size()),
    allFaces_
    (
        IOobject
        (
            "faces",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        faces
    ),
    faces_(allFaces_, allFaces_.size()),
    owner_
    (
        IOobject
        (
            "owner",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    neighbour_
    (
        IOobject
        (
            "neighbour",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0
    ),
    clearedPrimitives_(false),
    boundary_
    (
        IOobject
        (
            "boundary",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        *this,
        0
    ),
    bounds_(allPoints_, syncPar),
    directions_(Vector<label>::zero),
    pointZones_
    (
        IOobject
        (
            "pointZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    faceZones_
    (
        IOobject
        (
            "faceZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    cellZones_
    (
        IOobject
        (
            "cellZones",
            instance(),
            meshSubDir,
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    globalMeshDataPtr_(NULL),
    moving_(false),
    changing_(false),
    curMotionTimeIndex_(time().timeIndex()),
    oldAllPointsPtr_(NULL),
    oldPointsPtr_(NULL)
{
    // Check if the faces and cells are valid
    forAll (allFaces_, faceI)
    {
        const face& curFace = allFaces_[faceI];

        if (min(curFace) < 0 || max(curFace) > allPoints_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject& io,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const cellList& cells\n"
                ")\n"
            )   << "Face " << faceI << "contains vertex labels out of range: "
                << curFace << " Max point index = " << allPoints_.size()
                << abort(FatalError);
        }
    }

    // Check if the faces and cells are valid
    forAll (cells, cellI)
    {
        const cell& curCell = cells[cellI];

        if (min(curCell) < 0 || max(curCell) > allFaces_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh\n"
                "(\n"
                "    const IOobject& io,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const cellList& cells\n"
                ")\n"
            )   << "Cell " << cellI << "contains face labels out of range: "
                << curCell << " Max face index = " << allFaces_.size()
                << abort(FatalError);
        }
    }

    // Set the primitive mesh
    initMesh(const_cast<cellList&>(cells));
}


void Foam::polyMesh::resetPrimitives
(
    const label nUsedFaces,
    const pointField& points,
    const faceList& faces,
    const labelList& owner,
    const labelList& neighbour,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const bool validBoundary
)
{
    // Clear addressing. Keep geometric props for mapping.
    clearAddressing();

    // Take over new primitive data. Note extra optimization to prevent
    // assignment to self.
    if (&allPoints_ != &points)
    {
        allPoints_ = points;
        // Points will be reset in initMesh()
    }
    if (&allFaces_ != &faces)
    {
        allFaces_ = faces;
        // Faces will be reset in initMesh()
//         faces_.reset(allFaces_, owner.size());
    }
    if (&owner_ != &owner)
    {
        owner_ = owner;
    }
    if (&neighbour_ != &neighbour)
    {
        neighbour_ = neighbour;
    }

    // Reset patch sizes and starts
    forAll(boundary_, patchI)
    {
        boundary_[patchI] = polyPatch
        (
            boundary_[patchI].name(),
            patchSizes[patchI],
            patchStarts[patchI],
            patchI,
            boundary_
        );
    }


    // Flags the mesh files as being changed
    setInstance(time().timeName());

    // Check if the faces and cells are valid
    forAll (allFaces_, faceI)
    {
        const face& curFace = allFaces_[faceI];

        if (curFace.size() == 0)
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh::resetPrimitives\n"
                "(\n"
                "    const label nUsedFaces,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const labelList& owner,\n"
                "    const labelList& neighbour,\n"
                "    const labelList& patchSizes,\n"
                "    const labelList& patchStarts\n"
                ")\n"
            )   << "Face " << faceI << " contains no vertex labels"
                << abort(FatalError);
        }
        else if (min(curFace) < 0 || max(curFace) > allPoints_.size())
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh::resetPrimitives\n"
                "(\n"
                "    const label nUsedFaces,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const labelList& owner,\n"
                "    const labelList& neighbour,\n"
                "    const labelList& patchSizes,\n"
                "    const labelList& patchStarts\n"
                ")\n"
            )   << "Face " << faceI << " contains vertex labels out of range: "
                << curFace << " Max point index = " << allPoints_.size()
                << abort(FatalError);
        }
    }


    // Set the primitive mesh from the owner_, neighbour_. Works
    // out from patch end where the active faces stop.
    initMesh();

    // Recalculate bounds with all points.  HJ, 17/Oct/2008
    bounds_ = boundBox(allPoints_, validBoundary);

    if (validBoundary)
    {
        // Note that we assume that all the patches stay the same and are
        // correct etc. so we can already use the patches to do
        // processor-processor comms.

        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        // Warn if global empty mesh (constructs globalData!)
        if
        (
            globalData().nTotalPoints() == 0
         || globalData().nTotalCells() == 0
        )
        {
            FatalErrorIn
            (
                "polyMesh::polyMesh::resetPrimitives\n"
                "(\n"
                "    const label nUsedFaces,\n"
                "    const pointField& points,\n"
                "    const faceList& faces,\n"
                "    const labelList& owner,\n"
                "    const labelList& neighbour,\n"
                "    const labelList& patchSizes,\n"
                "    const labelList& patchStarts\n"
                ")\n"
            )
                << "no points or no cells in mesh" << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyMesh::~polyMesh()
{
    clearOut();
    resetMotion();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fileName& Foam::polyMesh::dbDir() const
{
    if (objectRegistry::dbDir() == defaultRegion)
    {
        return parent().dbDir();
    }
    else
    {
        return objectRegistry::dbDir();
    }
}


Foam::fileName Foam::polyMesh::meshDir() const
{
    return dbDir()/meshSubDir;
}


const Foam::fileName& Foam::polyMesh::pointsInstance() const
{
    return allPoints_.instance();
}


const Foam::fileName& Foam::polyMesh::facesInstance() const
{
    return allFaces_.instance();
}


const Foam::Vector<Foam::label>& Foam::polyMesh::directions() const
{
    if (directions_.x() == 0)
    {
        calcDirections();
    }

    return directions_;
}


Foam::label Foam::polyMesh::nGeometricD() const
{
    label nWedges = 0;

    forAll(boundary_, patchi)
    {
        if (isA<wedgePolyPatch>(boundary_[patchi]))
        {
            nWedges++;
        }
    }

    if (nWedges != 0 && nWedges != 2 && nWedges != 4)
    {
        FatalErrorIn("label polyMesh::nGeometricD() const")
            << "Number of wedge patches " << nWedges << " is incorrect, "
               "should be 0, 2 or 4"
            << exit(FatalError);
    }

    return nSolutionD() - nWedges/2;
}


Foam::label Foam::polyMesh::nSolutionD() const
{
    return cmptSum(directions() + Vector<label>::one)/2;
}


// Add boundary patches. Constructor helper
void Foam::polyMesh::addPatches
(
    const List<polyPatch*>& p,
    const bool validBoundary
)
{
    if (boundaryMesh().size() > 0)
    {
        FatalErrorIn
        (
            "void polyMesh::addPatches(const List<polyPatch*>&, const bool)"
        )   << "boundary already exists"
            << abort(FatalError);
    }

    boundary_.setSize(p.size());

    // Copy the patch pointers
    forAll (p, pI)
    {
        boundary_.set(pI, p[pI]);
    }

    // parallelData depends on the processorPatch ordering so force
    // recalculation. Problem: should really be done in removeBoundary but
    // there is some info in parallelData which might be interesting inbetween
    // removeBoundary and addPatches.
    deleteDemandDrivenData(globalMeshDataPtr_);

    if (validBoundary)
    {
        // Calculate topology for the patches (processor-processor comms etc.)
        boundary_.updateMesh();

        // Calculate the geometry for the patches (transformation tensors etc.)
        boundary_.calcGeometry();

        boundary_.checkDefinition();
    }
}


// Add mesh zones. Constructor helper
void Foam::polyMesh::addZones
(
    const List<pointZone*>& pz,
    const List<faceZone*>& fz,
    const List<cellZone*>& cz
)
{
    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        FatalErrorIn
        (
            "void addZones\n"
            "(\n"
            "    const List<pointZone*>& pz,\n"
            "    const List<faceZone*>& fz,\n"
            "    const List<cellZone*>& cz\n"
            ")"
        )   << "point, face or cell zone already exists"
            << abort(FatalError);
    }

    // Point zones
    if (pz.size())
    {
        pointZones_.setSize(pz.size());

        // Copy the zone pointers
        forAll (pz, pI)
        {
            pointZones_.set(pI, pz[pI]);
        }

        pointZones_.writeOpt() = IOobject::AUTO_WRITE;

        // Temporarily disable zone writing on creation
        // Auto-write should be sufficient.  HJ, 20/Aug/2009
//         pointZones_.write();
    }

    // Face zones
    if (fz.size())
    {
        faceZones_.setSize(fz.size());

        // Copy the zone pointers
        forAll (fz, fI)
        {
            faceZones_.set(fI, fz[fI]);
        }

        faceZones_.writeOpt() = IOobject::AUTO_WRITE;

        // Temporarily disable zone writing on creation
        // Auto-write should be sufficient.  HJ, 20/Aug/2009
//         faceZones_.write();
    }

    // Cell zones
    if (cz.size())
    {
        cellZones_.setSize(cz.size());

        // Copy the zone pointers
        forAll (cz, cI)
        {
            cellZones_.set(cI, cz[cI]);
        }

        cellZones_.writeOpt() = IOobject::AUTO_WRITE;

        // Temporarily disable zone writing on creation
        // Auto-write should be sufficient.  HJ, 20/Aug/2009
//         cellZones_.write();
    }
}


const Foam::pointField& Foam::polyMesh::allPoints() const
{
    if (clearedPrimitives_)
    {
        FatalErrorIn("const pointField& polyMesh::allPoints() const")
            << "allPoints deallocated"
            << abort(FatalError);
    }

    return allPoints_;
}


const Foam::pointField& Foam::polyMesh::points() const
{
    if (clearedPrimitives_)
    {
        FatalErrorIn("const pointField& polyMesh::points() const")
            << "points deallocated"
            << abort(FatalError);
    }

    return points_;
}


const Foam::faceList& Foam::polyMesh::allFaces() const
{
    if (clearedPrimitives_)
    {
        FatalErrorIn("const faceList& polyMesh::allFaces() const")
            << "allFaces deallocated"
            << abort(FatalError);
    }

    return allFaces_;
}


const Foam::faceList& Foam::polyMesh::faces() const
{
    if (clearedPrimitives_)
    {
        FatalErrorIn("const faceList& polyMesh::faces() const")
            << "allFaces deallocated"
            << abort(FatalError);
    }

    return faces_;
}


const Foam::labelList& Foam::polyMesh::faceOwner() const
{
    return owner_;
}


const Foam::labelList& Foam::polyMesh::faceNeighbour() const
{
    return neighbour_;
}


// Return old mesh motion points
const Foam::pointField& Foam::polyMesh::oldAllPoints() const
{
    if (!oldAllPointsPtr_)
    {
        if (debug)
        {
            WarningIn("const pointField& polyMesh::oldAllPoints() const")
                << "Old points not available.  Forcing storage of old points"
                << endl;
        }

        oldAllPointsPtr_ = new pointField(allPoints_);
        curMotionTimeIndex_ = time().timeIndex();
    }

    return *oldAllPointsPtr_;
}


// Return old mesh motion points
const Foam::pointField& Foam::polyMesh::oldPoints() const
{
    if (!oldPointsPtr_)
    {
        oldPointsPtr_ = new pointField::subField(oldAllPoints(), nPoints());
    }

    return *oldPointsPtr_;
}


Foam::tmp<Foam::scalarField> Foam::polyMesh::movePoints
(
    const pointField& newPoints
)
{
    if (debug)
    {
        Info<< "tmp<scalarField> polyMesh::movePoints(const pointField&) : "
            << " Moving points for time " << time().value()
            << " index " << time().timeIndex() << endl;
    }

    moving(true);

    // Pick up old points
    if (curMotionTimeIndex_ != time().timeIndex())
    {
        // Mesh motion in the new time step
        deleteDemandDrivenData(oldAllPointsPtr_);
        deleteDemandDrivenData(oldPointsPtr_);
        oldAllPointsPtr_ = new pointField(allPoints_);
        curMotionTimeIndex_ = time().timeIndex();
    }

    allPoints_ = newPoints;

    if (debug > 1)
    {
        // Check mesh motion
        if (primitiveMesh::checkMeshMotion(allPoints_, true))
        {
            Info<< "tmp<scalarField> polyMesh::movePoints"
                << "(const pointField&) : "
                << "Moving the mesh with given points will "
                << "invalidate the mesh." << nl
                << "Mesh motion should not be executed." << endl;
        }
    }

    allPoints_.writeOpt() = IOobject::AUTO_WRITE;
    allPoints_.instance() = time().timeName();

    points_.reset(allPoints_, nPoints());

    tmp<scalarField> sweptVols = primitiveMesh::movePoints
    (
        points_,
        oldAllPoints()
    );

    // Adjust parallel shared points
    if (globalMeshDataPtr_)
    {
        globalMeshDataPtr_->movePoints(allPoints_);
    }

    // Force recalculation of all geometric data with new points

    bounds_ = boundBox(allPoints_);
    boundary_.movePoints(allPoints_);

    pointZones_.movePoints(allPoints_);
    faceZones_.movePoints(allPoints_);
    cellZones_.movePoints(allPoints_);

    return sweptVols;
}


// Reset motion by deleting old points
void Foam::polyMesh::resetMotion() const
{
    curMotionTimeIndex_ = 0;
    deleteDemandDrivenData(oldAllPointsPtr_);
    deleteDemandDrivenData(oldPointsPtr_);
}


// Return parallel info
const Foam::globalMeshData& Foam::polyMesh::globalData() const
{
    if (!globalMeshDataPtr_)
    {
        if (debug)
        {
            Pout<< "polyMesh::globalData() const : "
                << "Constructing parallelData from processor topology"
                << endl;
        }
        // Construct globalMeshData using processorPatch information only.
        globalMeshDataPtr_ = new globalMeshData(*this);

        // Old method.  HJ, 6/Dec/2006

//         // Check for parallel boundaries
//         bool parBoundaries = false;

//         forAll (boundaryMesh(), patchI)
//         {
//             if
//             (
//                 typeid(boundaryMesh()[patchI])
//              == typeid(processorPolyPatch)
//             )
//             {
//                 parBoundaries = true;
//                 break;
//             }
//         }

//         if (parBoundaries)
//         {
//             // All is well - read the parallel data

//             globalDataPtr_ =
//                 new globalMeshData
//                 (
//                     IOobject
//                     (
//                         "globalData",
//                         time().findInstance(meshDir(), "globalData"),
//                         meshSubDir,
//                         *this,
//                         IOobject::MUST_READ,
//                         IOobject::NO_WRITE
//                     ),
//                     *this
//                 );
//         }
//         else
//         {
//             // The mesh has no parallel boundaries.  Create and hook a
//             // "non-parallel" parallel info
//             globalDataPtr_ =
//                 new globalMeshData
//                 (
//                     *this,
//                     false,
//                     false,  // cyclicParallel.  Remove when fixed
//                     nPoints(),
//                     nFaces(),
//                     nCells(),
//                     0,
//                     labelList(0),
//                     labelList(0),
//                     labelList(0)
//                 );
//         }
    }

    return *globalMeshDataPtr_;
}


// Remove all files and some subdirs (eg, sets)
void Foam::polyMesh::removeFiles(const fileName& instanceDir) const
{
    fileName meshFilesPath = db().path()/instanceDir/meshSubDir;

    rm(meshFilesPath/"points");
    rm(meshFilesPath/"faces");
    rm(meshFilesPath/"owner");
    rm(meshFilesPath/"neighbour");
    rm(meshFilesPath/"cells");
    rm(meshFilesPath/"boundary");
    rm(meshFilesPath/"pointZones");
    rm(meshFilesPath/"faceZones");
    rm(meshFilesPath/"cellZones");
    rm(meshFilesPath/"meshModifiers");
    rm(meshFilesPath/"parallelData");

    // remove subdirectories
    if (dir(meshFilesPath/"sets"))
    {
        rmDir(meshFilesPath/"sets");
    }
}

void Foam::polyMesh::removeFiles() const
{
    removeFiles(instance());
}


// ************************************************************************* //
