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

\*---------------------------------------------------------------------------*/

#include "faBoundaryMesh.H"
#include "faMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faBoundaryMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
faBoundaryMesh::faBoundaryMesh
(
    const IOobject& io,
    const faMesh& mesh
)
:
    faPatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    if (readOpt() == IOobject::MUST_READ)
    {
        faPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchI)
        {
            patches.set
            (
                patchI,
                faPatch::New
                (
                    patchEntries[patchI].keyword(),
                    patchEntries[patchI].dict(),
                    patchI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "faBoundaryMesh::polyBoundaryMesh"
            "(const IOobject&, const faMesh&)"
        );

        close();
    }
}


// Construct given size. Patches will be set later
faBoundaryMesh::faBoundaryMesh
(
    const IOobject& io,
    const faMesh& pm,
    const label size
)
:
    faPatchList(size),
    regIOobject(io),
    mesh_(pm)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate the geometry for the patches (transformation tensors etc.)
void faBoundaryMesh::calcGeometry()
{
    forAll(*this, patchi)
    {
        operator[](patchi).initGeometry();
    }

    forAll(*this, patchi)
    {
        operator[](patchi).calcGeometry();
    }
}


// Return the mesh reference
const faMesh& faBoundaryMesh::mesh() const
{
    return mesh_;
}


lduInterfacePtrsList faBoundaryMesh::interfaces() const
{
    lduInterfacePtrsList interfaces(size());

    forAll (interfaces, patchi)
    {
        if (isA<lduInterface>(this->operator[](patchi)))
        {
            interfaces.set
            (
                patchi,
               &refCast<const lduInterface>(this->operator[](patchi))
            );
        }
    }

    return interfaces;
}


// Return a list of patch types
wordList faBoundaryMesh::types() const
{
    const faPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].type();
    }

    return t;
}


// Return a list of patch names
wordList faBoundaryMesh::names() const
{
    const faPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].name();
    }

    return t;
}


label faBoundaryMesh::findPatchID(const word& patchName) const
{
    const faPatchList& patches = *this;

    forAll (patches, patchI)
    {
        if (patches[patchI].name() == patchName)
        {
            return patchI;
        }
    }

    // Patch not found
    return -1;
}


// Return patch index for a given edge label
label faBoundaryMesh::whichPatch(const label edgeIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (edgeIndex >= mesh().nEdges())
    {
        FatalErrorIn
        (
            "faBoundaryMesh::whichPatch(const label edgeIndex) const"
        )   << "given label greater than the number of edges"
            << abort(FatalError);
    }

    if (edgeIndex < mesh().nInternalEdges())
    {
        return -1;
    }

    forAll (*this, patchI)
    {
        const faPatch& bp = operator[](patchI);

        if
        (
            edgeIndex >= bp.start()
         && edgeIndex < bp.start() + bp.size()
        )
        {
            return patchI;
        }
    }

    // If not in any of above, it's trouble!
    FatalErrorIn
    (
        "label faBoundaryMesh::whichPatch(const label edgeIndex) const"
    )   << "error in patch search algorithm"
        << abort(FatalError);

    return -1;
}


bool faBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalEdges();
    const faBoundaryMesh& bm = *this;

    bool boundaryError = false;

    forAll (bm, patchI)
    {
        if (bm[patchI].start() != nextPatchStart)
        {
            boundaryError = true;

            Info
                << "bool faBoundaryMesh::checkDefinition("
                << "const bool report) const : "
                << "Problem with boundary patch " << patchI
                << ".\nThe patch should start on face no " << nextPatchStart
                << " and the boundary file specifies " << bm[patchI].start()
                << "." << nl << endl;
        }

        nextPatchStart += bm[patchI].faPatch::size();
    }

    if (boundaryError)
    {
        SeriousErrorIn
        (
            "bool faBoundaryMesh::checkDefinition("
            "const bool report) const"
        )   << "This mesh is not valid: boundary definition is in error."
            << endl;
    }
    else
    {
        if (debug || report)
        {
            Info << "Boundary definition OK." << endl;
        }
    }

    return boundaryError;
}


// Correct faBoundaryMesh after moving points
void faBoundaryMesh::movePoints(const pointField& p)
{
    faPatchList& patches = *this;

    forAll (patches, patchI)
    {
        patches[patchI].initMovePoints(p);
    }

    forAll (patches, patchI)
    {
        patches[patchI].movePoints(p);
    }
}


void faBoundaryMesh::updateMesh()
{
    faPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initUpdateMesh();
    }

    forAll(patches, patchi)
    {
        patches[patchi].updateMesh();
    }
}


// writeData member function required by regIOobject
bool faBoundaryMesh::writeData(Ostream& os) const
{
    const faPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches, patchi)
    {
        os  << indent << patches[patchi].name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << patches[patchi] << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_LIST;

    // Check state of IOstream
    os.check("polyBoundaryMesh::writeData(Ostream& os) const");

    return os.good();
}


Ostream& operator<<(Ostream& os, const faBoundaryMesh& bm)
{
    bm.writeData(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
