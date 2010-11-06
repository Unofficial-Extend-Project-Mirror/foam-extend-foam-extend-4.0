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

Description

    Mesh needed to do edge-based addressing.

\*---------------------------------------------------------------------------*/

#include "eMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(eMesh, 0);

word eMesh::meshSubDir = "eMesh";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void eMesh::clearGeom() const
{
    if (debug)
    {
        Info<< "void eMesh::clearGeom() const : "
            << "Clearing geometry" << endl;
    }
}


void eMesh::clearAddressing() const
{
    if (debug)
    {
        Info<< "void eMesh::clearAddressing() const : "
            << "Clearing addressing" << endl;
    }

    edges_.clear();

    deleteDemandDrivenData(pePtr_);
    deleteDemandDrivenData(epPtr_);
    deleteDemandDrivenData(fePtr_);
    deleteDemandDrivenData(efPtr_);
}


// Helper function to isolate points on triangular faces
label eMesh::findIsolatedPoint(const face& f, const edge& e) const
{
    // Check the first point
    if ( f[0] != e.start() && f[0] != e.end() )
    {
        return f[0];
    }

    // Check the second point
    if ( f[1] != e.start() && f[1] != e.end() )
    {
        return f[1];
    }

    // Check the third point
    if ( f[2] != e.start() && f[2] != e.end() )
    {
        return f[2];
    }

    return -1;
}


//- Helper function to determine the orientation of a triangular face
label eMesh::edgeDirection(const face& f, const edge& e) const
{
    if
    (
        (f[0] == e.start() && f[1] == e.end())
     || (f[1] == e.start() && f[2] == e.end())
     || (f[2] == e.start() && f[0] == e.end())
    )
    {
        // Return counter-clockwise
        return 1;
    }
    else
    {
        // Return clockwise
        return -1;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

eMesh::eMesh(const polyMesh& pMesh, const word& subDir)
:
    objectRegistry(pMesh.time()),
    mesh_(pMesh),
    edges_
    (
        IOobject
        (
            "edges",
            mesh_.facesInstance(),
            subDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    boundary_
    (
        IOobject
        (
            "edgeBoundary",
            mesh_.facesInstance(),
            subDir,
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        *this
    ),
    pePtr_(NULL),
    epPtr_(NULL),
    fePtr_(NULL),
    efPtr_(NULL)
{
    if (debug)
    {
        Info << "eMesh::eMesh(const polyMesh&, const word&) : "
             << "Creating eMesh from polyMesh"
             << endl;
    }

    // Re-initialize / override meshSubDir
    meshSubDir = subDir;

    // Try to read from disk.
    if (edges_.headerOk() && boundary_.headerOk())
    {
        // Set sizes
        nEdges_ = edges_.size();
        nInternalEdges_ = boundary_[0].start();
    }
    else
    {
        // Could not read ordered edges, so calculate it instead.
        calcOrderedEdgeList();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

eMesh::~eMesh()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

fileName eMesh::meshDir() const
{
    return mesh_.dbDir();
}


fileName eMesh::meshSubDirectory() const
{
    return mesh_.dbDir()/meshSubDir;
}


const Time& eMesh::time() const
{
    return mesh_.time();
}


label eMesh::nEdges() const
{
    return nEdges_;
}


label eMesh::nInternalEdges() const
{
    return nInternalEdges_;
}


const edgeList& eMesh::edges() const
{
    return edges_;
}


void eMesh::addEdgePatches(const List<ePatch*>& p)
{
    if (debug)
    {
        Info << "void eMesh::addEdgePatches(const List<ePatch*>& p) : "
             << "Adding patches to eMesh" << endl;
    }

    if (boundary().size() > 0)
    {
        FatalErrorIn
        (
            "void eMesh::addEdgePatches(const List<ePatch*>& p)"
        )
            << "Boundary already exists"
            << abort(FatalError);
    }

    boundary_.setSize(p.size());

    forAll(p, patchI)
    {
        boundary_.set(patchI, p[patchI]);
    }
}


const objectRegistry& eMesh::db() const
{
    return mesh_.db();
}


const eBoundaryMesh& eMesh::boundary() const
{
    return boundary_;
}


void eMesh::resetPrimitives
(
    edgeList& edges,
    labelListList& faceEdges,
    labelListList& edgeFaces,
    const labelList& patchSizes,
    const labelList& patchStarts,
    const bool reUse,
    const bool storePrimitives
)
{
    // Clear out geometry and addressing
    clearOut();

    // Initialize pointers for storage
    if (storePrimitives)
    {
        fePtr_ =
        (
            new labelListIOList
            (
                IOobject
                (
                    "faceEdges",
                    mesh_.facesInstance(),
                    meshSubDir,
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        efPtr_ =
        (
            new labelListIOList
            (
                IOobject
                (
                    "edgeFaces",
                    mesh_.facesInstance(),
                    meshSubDir,
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        if (reUse)
        {
            edges_.transfer(edges);
            fePtr_->transfer(faceEdges);
            efPtr_->transfer(edgeFaces);
        }
        else
        {
            edges_ = edges;
            fePtr_->operator=(faceEdges);
            efPtr_->operator=(edgeFaces);
        }
    }

    // Reset patch sizes and starts
    forAll(boundary_, patchI)
    {
        boundary_[patchI] = ePatch
        (
            boundary_[patchI].name(),
            patchSizes[patchI],
            patchStarts[patchI],
            patchI,
            boundary_
        );
    }

    // Reset size information
    nEdges_ = edges.size();
    nInternalEdges_ = boundary_[0].start();

    // Set mesh files as changed
    setInstance(time().timeName());
}


//- Clear demand-driven data
void eMesh::clearOut() const
{
    clearGeom();
    clearAddressing();
}


//- Set the instance for mesh files
void eMesh::setInstance(const fileName& inst)
{
    if (debug)
    {
        Info<< "void eMesh::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    if (edges_.size())
    {
        edges_.writeOpt() = IOobject::AUTO_WRITE;
        edges_.instance() = inst;
    }

    if (efPtr_)
    {
        efPtr_->writeOpt() = IOobject::AUTO_WRITE;
        efPtr_->instance() = inst;
    }

    if (fePtr_)
    {
        fePtr_->writeOpt() = IOobject::AUTO_WRITE;
        fePtr_->instance() = inst;
    }

    // Write out boundary information only
    // if all others are being written out
    if (edges_.size() && efPtr_ && fePtr_)
    {
        boundary_.writeOpt() = IOobject::AUTO_WRITE;
        boundary_.instance() = inst;
    }
}


const labelListList& eMesh::pointEdges() const
{
    if (!pePtr_)
    {
        calcPointEdges();
    }

    return *pePtr_;
}


const labelListList& eMesh::edgePoints() const
{
    if (!epPtr_)
    {
        calcEdgePoints();
    }

    return *epPtr_;
}


const labelListList& eMesh::faceEdges() const
{
    if (!fePtr_)
    {
        calcFaceEdges();
    }

    return *fePtr_;
}


const labelListList& eMesh::edgeFaces() const
{
    if (!efPtr_)
    {
        calcEdgeFaces();
    }

    return *efPtr_;
}


bool eMesh::write() const
{
    if (edges_.size())
    {
        edges_.write();
    }

    if (efPtr_)
    {
        efPtr_->write();
    }

    if (fePtr_)
    {
        fePtr_->write();
    }

    if (edges_.size() && efPtr_ && fePtr_)
    {
        boundary_.write();
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool eMesh::operator!=(const eMesh& m) const
{
    return &m != this;
}

bool eMesh::operator==(const eMesh& m) const
{
    return &m == this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
