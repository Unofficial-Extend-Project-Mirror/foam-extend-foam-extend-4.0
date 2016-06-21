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

Class
    Foam::polyBoundaryMesh

Description
    Foam::polyBoundaryMesh

SourceFiles
    polyBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef polyBoundaryMesh_H
#define polyBoundaryMesh_H

#include "polyPatchList.H"
#include "regIOobject.H"
#include "labelPairList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

// Forward declaration of friend functions and operators

class polyBoundaryMesh;

Ostream& operator<<(Ostream&, const polyBoundaryMesh&);


/*---------------------------------------------------------------------------*\
                      Class polyBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class polyBoundaryMesh
:
    public polyPatchList,
    public regIOobject
{
    // private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Edges of neighbouring patches
        mutable labelPairListList* neighbourEdgesPtr_;


    // Private Member Functions

        //- Create identity map
        static labelList ident(const label len);

        //- Calculate the geometry for the patches,
        //  eg transformation tensors etc.
        //  Enhanced for parallel addressing calculations
        //  HJ, 11/Jun/2011
        void calcGeometry();

        //- Disallow construct as copy
        polyBoundaryMesh(const polyBoundaryMesh&);

        //- Disallow assignment
        void operator=(const polyBoundaryMesh&);


public:

    //- Declare friendship with polyMesh
    friend class polyMesh;


    //- Runtime type information
    TypeName("polyBoundaryMesh");


    // Constructors

        //- Read constructor given IOobject and a polyMesh reference
        //  Note point pointers are unset, only used in copying meshes
        polyBoundaryMesh
        (
            const IOobject&,
            const polyMesh&
        );

        //- Construct given size
        polyBoundaryMesh
        (
            const IOobject&,
            const polyMesh&,
            const label size
        );


    // Destructor

        ~polyBoundaryMesh();

        //- Clear geometry at this level and at patches
        void clearGeom();

        //- Clear addressing at this level and at patches
        void clearAddressing();


    // Member functions

        //- Return the mesh reference
        const polyMesh& mesh() const
        {
            return mesh_;
        }

        //- Per patch the edges on the neighbouring patch
        //  For every external edge the neighbouring patch and
        //  neighbouring (external) patch edge label.  Note that edge
        //  indices are offset by nInternalEdges to keep it as much as
        //  possible consistent with coupled patch addressing
        //  (where coupling is by local patch face index).
        //  Only valid for singly connected polyBoundaryMesh and not parallel
        const labelPairListList& neighbourEdges() const;

        //- Return a list of patch names
        wordList names() const;

        //- Return a list of patch types
        wordList types() const;

        //- Return a list of physical types
        wordList physicalTypes() const;

        //- Find patch index given a name
        label findPatchID(const word& patchName) const;

        //- Return patch index for a given face label
        label whichPatch(const label faceIndex) const;

        //- Return the set of patch IDs corresponding to the given list
        //  of names.  Wild cards are expanded.
        labelHashSet patchSet(const wordList&) const;

        //- Check whether all procs have all patches and in same order. Return
        //  true if in error.
        bool checkParallelSync(const bool report = false) const;

        //- Check boundary definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Correct polyBoundaryMesh after moving points
        void movePoints(const pointField&);

        //- Correct after topology update  Note: Boundary has been rebuilt
        //  and no mapping is required.  HJ, 1/Mar/2009
        void updateMesh();

        //- Reorders patches. Ordering does not have to be done in
        //  ascending or descending order. Reordering has to be unique.
        //  (is shuffle) Calls updateMesh() after reordering to recalculate
        //  data.
        void reorder(const UList<label>&);

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Write using given format, version and form uncompression
        bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const;


    // Member Operators

        //- Return const and non-const reference to polyPatch by index.
        using polyPatchList::operator[];

        //- Return const reference to polyPatch by name.
        const polyPatch& operator[](const word&) const;

        //- Return reference to polyPatch by name.
        polyPatch& operator[](const word&);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const polyBoundaryMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
