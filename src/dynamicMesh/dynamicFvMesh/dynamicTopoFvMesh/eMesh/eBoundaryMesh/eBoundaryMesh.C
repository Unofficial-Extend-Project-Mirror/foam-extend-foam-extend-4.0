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

#include "eBoundaryMesh.H"
#include "eMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(eBoundaryMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
eBoundaryMesh::eBoundaryMesh
(
    const IOobject& io,
    const eMesh& mesh
)
:
    ePatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    if
    (
        (readOpt() == IOobject::MUST_READ) ||
        (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readFromInputStream();
    }
}


// Construct given size. Patches will be set later
eBoundaryMesh::eBoundaryMesh
(
    const IOobject& io,
    const eMesh& em,
    const label size
)
:
    ePatchList(size),
    regIOobject(io),
    mesh_(em)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Read from input stream
void eBoundaryMesh::readFromInputStream()
{
    ePatchList& patches = *this;

    // Read polyPatchList
    Istream& is = readStream(typeName);

    PtrList<entry> patchEntries(is);
    patches.setSize(patchEntries.size());

    forAll(patches, patchI)
    {
        patches.set
        (
            patchI,
            ePatch::New
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
        "eBoundaryMesh::eBoundaryMesh"
        "(const IOobject&, const faMesh&)"
    );

    close();
}

// Return the mesh reference
const eMesh& eBoundaryMesh::mesh() const
{
    return mesh_;
}


// Return a list of patch types
wordList eBoundaryMesh::types() const
{
    const ePatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].type();
    }

    return t;
}


// Return a list of patch names
wordList eBoundaryMesh::names() const
{
    const ePatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].name();
    }

    return t;
}


// Return a list of patch sizes
labelList eBoundaryMesh::patchSizes() const
{
    const ePatchList& patches = *this;

    labelList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].size();
    }

    return t;
}


// Return a list of patch starts
labelList eBoundaryMesh::patchStarts() const
{
    const ePatchList& patches = *this;

    labelList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].start();
    }

    return t;
}


//- Return patch index for a given edge label
label eBoundaryMesh::whichPatch(const label edgeIndex) const
{
    // Find out which patch the current edge belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (edgeIndex >= mesh().nEdges())
    {
        FatalErrorIn
        (
            "eBoundaryMesh::whichPatch(const label edgeIndex) const"
        )   << "given label greater than the number of geometric edges"
            << abort(FatalError);
    }

    if (edgeIndex < mesh().nInternalEdges())
    {
        return -1;
    }

    forAll (*this, patchI)
    {
        const ePatch& bp = operator[](patchI);

        if
        (
            edgeIndex >= bp.start()
         && edgeIndex < bp.start() + bp.size()
        )
        {
            return patchI;
        }
    }

    // If not in any of above, it is trouble!
    FatalErrorIn
    (
        "label eBoundaryMesh::whichPatch(const label edgeIndex) const"
    )   << "Cannot find edge " << edgeIndex << " in any of the patches "
        << names() << nl
        << "It seems your patches are not consistent with the mesh :"
        << " internalEdges:" << mesh().nInternalEdges()
        << "  total number of edges:" << mesh().nEdges()
        << abort(FatalError);

    return -1;
}


label eBoundaryMesh::findPatchID(const word& patchName) const
{
    const ePatchList& patches = *this;

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


// writeData member function required by regIOobject
bool eBoundaryMesh::writeData(Ostream& os) const
{
    const ePatchList& patches = *this;

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
    os.check("eBoundaryMesh::writeData(Ostream& os) const");

    return os.good();
}


Ostream& operator<<(Ostream& os, const eBoundaryMesh& bm)
{
    bm.writeData(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
