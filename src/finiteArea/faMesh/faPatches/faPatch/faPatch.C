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

#include "faPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faPatch, 0);
    defineRunTimeSelectionTable(faPatch, dictionary);
    addToRunTimeSelectionTable(faPatch, faPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faPatch::clearOut()
{
    deleteDemandDrivenData(edgeFacesPtr_);
    deleteDemandDrivenData(pointLabelsPtr_);
    deleteDemandDrivenData(pointEdgesPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faPatch::faPatch
(
    const word& name,
    const labelList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm,
    const label ngbPolyPatchIndex
)
:
    labelList(edgeLabels),
    patchIdentifier(name, index),
    ngbPolyPatchIndex_(ngbPolyPatchIndex),
    boundaryMesh_(bm),
    edgeFacesPtr_(NULL),
    pointLabelsPtr_(NULL),
    pointEdgesPtr_(NULL)
{}


// Construct from dictionary
Foam::faPatch::faPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
:
    labelList(dict.lookup("edgeLabels")),
    patchIdentifier(name, dict, index),
    ngbPolyPatchIndex_(readInt(dict.lookup("ngbPolyPatchIndex"))),
    boundaryMesh_(bm),
    edgeFacesPtr_(NULL),
    pointLabelsPtr_(NULL),
    pointEdgesPtr_(NULL)
{}

Foam::faPatch::faPatch(const faPatch& p, const faBoundaryMesh& bm)
:
    labelList(p),
    patchIdentifier(p, p.index()),
    ngbPolyPatchIndex_(p.ngbPolyPatchIndex_),
    boundaryMesh_(bm),
    edgeFacesPtr_(NULL),
    pointLabelsPtr_(NULL),
    pointEdgesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faPatch::~faPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::faPatch::ngbPolyPatchIndex() const
{
    return ngbPolyPatchIndex_;
}

const Foam::faBoundaryMesh& Foam::faPatch::boundaryMesh() const
{
    return boundaryMesh_;
}


Foam::label Foam::faPatch::start() const
{
    return boundaryMesh().mesh().patchStarts()[index()];
}


const Foam::labelList& Foam::faPatch::pointLabels() const
{
    if (!pointLabelsPtr_)
    {
        calcPointLabels();
    }

    return *pointLabelsPtr_;
}


void Foam::faPatch::calcPointLabels() const
{
    SLList<label> labels;

    UList<edge> edges =
        patchSlice(boundaryMesh().mesh().edges());

    forAll(edges, edgeI)
    {
        bool existStart = false;
        bool existEnd = false;

        for
        (
            SLList<label>::iterator iter = labels.begin();
            iter != labels.end();
            ++iter
        )
        {
            if(*iter == edges[edgeI].start())
            {
                existStart = true;
            }

            if(*iter == edges[edgeI].end())
            {
                existEnd = true;
            }
        }

        if(!existStart)
        {
            labels.append(edges[edgeI].start());
        }

        if(!existEnd)
        {
            labels.append(edges[edgeI].end());
        }
    }

    pointLabelsPtr_ = new labelList(labels);
}


void Foam::faPatch::calcPointEdges() const
{
    labelList points = pointLabels();

    const edgeList::subList e =
        patchSlice(boundaryMesh().mesh().edges());

    // set up storage for pointEdges
    List<SLList<label> > pointEdgs(points.size());

    forAll (e, edgeI)
    {
        const edge& curPoints = e[edgeI];

        forAll (curPoints, pointI)
        {
            label localPointIndex =
                findIndex(points, curPoints[pointI]);

            pointEdgs[localPointIndex].append(edgeI);
        }
    }

    // sort out the list
    pointEdgesPtr_ = new labelListList(pointEdgs.size());
    labelListList& pEdges = *pointEdgesPtr_;

    forAll (pointEdgs, pointI)
    {
        pEdges[pointI].setSize(pointEdgs[pointI].size());

        label i = 0;
        for
        (
            SLList<label>::iterator curEdgesIter = pointEdgs[pointI].begin();
            curEdgesIter != pointEdgs[pointI].end();
            ++curEdgesIter, ++i
        )
        {
            pEdges[pointI][i] = curEdgesIter();
        }
    }
}


const Foam::labelListList& Foam::faPatch::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        calcPointEdges();
    }

    return *pointEdgesPtr_;
}


Foam::labelList Foam::faPatch::ngbPolyPatchFaces() const
{
    labelList ngbFaces;

    if(ngbPolyPatchIndex() == -1)
    {
        return ngbFaces;
    }

    ngbFaces.setSize(faPatch::size());

    const faMesh& aMesh = boundaryMesh().mesh();
    const polyMesh& pMesh = aMesh();
    const indirectPrimitivePatch& patch = aMesh.patch();

    const labelListList& edgeFaces = pMesh.edgeFaces();

    labelList faceCells (patch.size(), -1);

    forAll (faceCells, faceI)
    {
        label faceID = aMesh.faceLabels()[faceI];

        faceCells[faceI] = pMesh.faceOwner()[faceID];
    }

    labelList meshEdges =
        patch.meshEdges
        (
            pMesh.edges(),
            pMesh.cellEdges(),
            faceCells
        );

    forAll(ngbFaces, edgeI)
    {
        ngbFaces[edgeI] = -1;

        label curEdge = (*this)[edgeI];

        label curPMeshEdge = meshEdges[curEdge];

        forAll(edgeFaces[curPMeshEdge], faceI)
        {
            label curFace = edgeFaces[curPMeshEdge][faceI];

            label curPatchID = pMesh.boundaryMesh().whichPatch(curFace);

            if (curPatchID == ngbPolyPatchIndex())
            {
                ngbFaces[edgeI] = curFace;
            }
        }

        if(ngbFaces[edgeI] == -1)
        {
            Info<< "faPatch::edgeNgbPolyPatchFaces(): "
                << "Problem with determination of edge ngb faces!" << endl;
        }
    }

    return ngbFaces;
}


Foam::tmp<Foam::vectorField> Foam::faPatch::ngbPolyPatchFaceNormals() const
{
    tmp<vectorField> tfN(new vectorField());
    vectorField& fN = tfN();

    if (ngbPolyPatchIndex() == -1)
    {
        return tfN;
    }

    fN.setSize(faPatch::size());

    labelList ngbFaces = ngbPolyPatchFaces();

    const polyMesh& pMesh = boundaryMesh().mesh()();

    const faceList& faces = pMesh.faces();
    const pointField& points = pMesh.points();

    forAll(fN, faceI)
    {
        fN[faceI] = faces[ngbFaces[faceI]].normal(points)
            /faces[ngbFaces[faceI]].mag(points);
    }

    return tfN;
}


Foam::tmp<Foam::vectorField> Foam::faPatch::ngbPolyPatchPointNormals() const
{
    if (ngbPolyPatchIndex() == -1)
    {
        return tmp<vectorField>(new vectorField());
    }

    labelListList pntEdges = pointEdges();

    tmp<vectorField> tpN(new vectorField(pntEdges.size(), vector::zero));
    vectorField& pN = tpN();

    vectorField faceNormals = ngbPolyPatchFaceNormals();

    forAll(pN, pointI)
    {
        forAll(pntEdges[pointI], edgeI)
        {
            pN[pointI] += faceNormals[pntEdges[pointI][edgeI]];
        }
    }

    pN /= mag(pN);

    return tpN;
}


const Foam::unallocLabelList& Foam::faPatch::edgeFaces() const
{
    if (!edgeFacesPtr_)
    {
        edgeFacesPtr_ = new labelList::subList
        (
            patchSlice(boundaryMesh().mesh().edgeOwner())
        );
    }

    return *edgeFacesPtr_;
}


// Return the patch edge centres
const Foam::vectorField& Foam::faPatch::edgeCentres() const
{
    return boundaryMesh().mesh().edgeCentres().boundaryField()[index()];
}


// Return the patch edges length vectors
const Foam::vectorField& Foam::faPatch::edgeLengths() const
{
    return boundaryMesh().mesh().Le().boundaryField()[index()];
}


// Return the patch edge length magnitudes
const Foam::scalarField& Foam::faPatch::magEdgeLengths() const
{
    return boundaryMesh().mesh().magLe().boundaryField()[index()];
}


// Return the patch edge unit normals
Foam::tmp<Foam::vectorField> Foam::faPatch::edgeNormals() const
{
    tmp<vectorField> eN(new vectorField(size()));

    eN() = edgeLengths()/magEdgeLengths();

    return eN;
}


// Return the patch edge neighbour face centres
Foam::tmp<Foam::vectorField> Foam::faPatch::edgeFaceCentres() const
{
    tmp<vectorField> tfc(new vectorField(size()));
    vectorField& fc = tfc();

    // get reference to global face centres
    const vectorField& gfc =
        boundaryMesh().mesh().areaCentres().internalField();

    const unallocLabelList& faceLabels = edgeFaces();

    forAll (faceLabels, edgeI)
    {
        fc[edgeI] = gfc[faceLabels[edgeI]];
    }

    return tfc;
}


// Return cell-centre to face-centre vector
Foam::tmp<Foam::vectorField> Foam::faPatch::delta() const
{
    return edgeCentres() - edgeFaceCentres();
}


// Make delta coefficients as patch face - neighbour cell distances
void Foam::faPatch::makeDeltaCoeffs(scalarField& dc) const
{
    dc = 1.0/(edgeNormals() & delta());
}


// Return delta coefficients
const Foam::scalarField& Foam::faPatch::deltaCoeffs() const
{
    return boundaryMesh().mesh().deltaCoeffs().boundaryField()[index()];
}


void Foam::faPatch::makeWeights(scalarField& w) const
{
    w = 1.0;
}


const Foam::scalarField& Foam::faPatch::weights() const
{
    return boundaryMesh().mesh().weights().boundaryField()[index()];
}


void Foam::faPatch::movePoints(const pointField& points)
{}


void Foam::faPatch::resetEdges(const labelList& newEdges)
{
    Info<< "Resetting patch edges" << endl;
    labelList::operator=(newEdges);

    clearOut();
}


void Foam::faPatch::write(Ostream& os) const
{
    os.writeKeyword("type") << type() << token::END_STATEMENT << nl;
    patchIdentifier::write(os);

    const labelList& edgeLabels = *this;
    edgeLabels.writeEntry("edgeLabels", os);
    os.writeKeyword("ngbPolyPatchIndex") << ngbPolyPatchIndex_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faPatch& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const faPatch& p)");
    return os;
}


// ************************************************************************* //
