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

Class
    Foam::meshCutSurface

Description
    Description of surface resulting from cuts through mesh. Region on
    labelledTris is the cell the triangle is resulting from.

SourceFiles
    meshCutSurface.C

\*---------------------------------------------------------------------------*/

#ifndef meshCutSurface_H
#define meshCutSurface_H

#include "labelList.H"
#include "HashTable.H"
#include "pointField.H"
#include "triSurface.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;
class cellAddressing;
class faceDecompCuts;
class cellDecompCuts;

/*---------------------------------------------------------------------------*\
                           Class meshCutSurface Declaration
\*---------------------------------------------------------------------------*/

class meshCutSurface
:
    public triSurface
{

    // Static data members

        //- Symbolic names of tet vertices
        static label FP0;
        static label FC;
        static label FP1;
        static label CC;

        //- Symbolic names of tet edges
        static label FP0_FC;
        static label FP0_FP1;
        static label FP0_CC;
        static label FC_CC;
        static label FC_FP1;
        static label FP1_CC;


    // Private Static Functions

        //- Find key in hashtable. Returns value or -1 if not found
        template<class T, class Hash>
        static label find
        (
            const HashTable<label, T, Hash>& table,
            const T& key
        );

        //- Find edge in subset edgeLabels of edges.
        static label findEdge
        (
            const edgeList& edges,
            const labelList& edgeLabels,
            const edge& e
        );

        //- Find the other (i.e. not cutEdgeI) cut-edge on faceI of cell model
        static label findCutEdge
        (
            const cellAddressing& model,
            const labelList& tetEdgeCuts,
            const label faceI,
            const label cutEdgeI
        );


    // Private Member Functions

        //- Remove duplicate faces (resulting from 3 meshVerts cut - is seen
        //  from both sides)
        static triSurface removeDuplicates(const triSurface& surf);


public:

        //- Generic tet cutting routine.
        //  Generate triangles given the cut vertices and edges of the tet.
        //  (Triangles can be incorrectly oriented)
        void cutTet
        (
            const cellAddressing& model,
            const labelList& tetVertCuts,
            const labelList& tetEdgeCuts,
            const label cellI,
            List<labelledTri>& tris,
            label& triI
        );

        //- Tet cutting when cuts are only through edges.
        void cutTetThroughEdges
        (
            const cellAddressing& model,
            const labelList& tetEdgeCuts,
            const label nEdgeCuts,
            const label cellI,
            List<labelledTri>& tris,
            label& triI
        );

        //- Tet cutting when cuts are through corner vertices (and possibly
        //  edges)
        void cutTetThroughVerts
        (
            const cellAddressing& model,
            const labelList& tetVertCuts,
            const label nVertCuts,
            const labelList& tetEdgeCuts,
            const label nEdgeCuts,
            const label cellI,
            List<labelledTri>& tris,
            label& triI
        );


    // Static functions

        //- Interpolate onto points given edge weights, cell centre values,
        //  face centre values and vertex values
        template <class T>
        static tmp<Field<T> > interpolate
        (
            const faceDecompCuts& cuts,
            const Field<T>& vField,
            const Field<T>& fField,
            const Field<T>& pField
        );

        //- Interpolate onto points given edge weights, cell centre values,
        //  and vertex values
        template <class T>
        static tmp<Field<T> > interpolate
        (
            const cellDecompCuts& cuts,
            const Field<T>& vField,
            const Field<T>& pField
        );


    // Constructors

        //- Construct from cuts on faceDecomp tet decomposition
        //  (introduces cellcentres and facecentres)
        meshCutSurface(const faceDecompCuts& cuts);

        //- Construct from cuts on cellDecomp tet decomposition
        //  (introduces cellcentres but not facecentres)
        meshCutSurface(const cellDecompCuts& cuts);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "meshCutSurfaceInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
