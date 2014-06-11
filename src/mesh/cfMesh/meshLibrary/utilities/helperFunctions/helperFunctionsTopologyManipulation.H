/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    helperFunctions

Description
    Topology modifiers and queries useful for mesh generation

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef helperFunctionsTopologyManipulation_H
#define helperFunctionsTopologyManipulation_H

#include "DynList.H"
#include "face.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class boolList;
class cellList;
class labelHashSet;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Namespace help functions Declaration
\*---------------------------------------------------------------------------*/

namespace help
{
    //- check if the faces are equal
    template<class faceType1, class faceType2>
    bool areFacesEqual(const faceType1& f1, const faceType2& f2);

    //- local position of element in a list
    template<class T, class ListType>
    label positionInList(const T& elmt, const ListType& l);

    //- reverse the face
    template<class faceType>
    faceType reverseFace(const faceType& f);

    //- returns a merged face
    template<class faceType1, class faceType2>
    inline face mergeTwoFaces
    (
        const faceType1& f1,
        const faceType2& f2
    );

    //- remove edges until faces share only one edge
    inline edgeList modifyFacesToShareOneEdge(face& f1, face& f2);

    //- create a face from the removed part
    inline face createFaceFromRemovedPart(const face& fOrig, const face& fCut);

    //- remove edges from face
    inline face removeEdgesFromFace
    (
        const face& fOrig,
        const DynList<edge>& removeEdges
    );

    //- find open edges for a set of faces forming a cell
    inline void findOpenEdges
    (
        const faceList& cellFaces,
        DynList<edge>& openEdges
    );

    //- check if two faces share an edge
    template<class faceType1, class faceType2>
    inline bool shareAnEdge(const faceType1& f1, const faceType2& f2);

    //- return the edge shared by the faces
    template<class faceType1, class faceType2>
    inline edge sharedEdge(const faceType1& f1, const faceType2& f2);

    //- return the position of edge in the face, -1 otherwise
    template<class faceType>
    inline label positionOfEdgeInFace(const edge& e, const faceType& f);

    //- check if two faces share a vertex
    template<class faceType1, class faceType2>
    inline bool shareAVertex(const faceType1& f1, const faceType2& f2);

    //- shared vertex of two faces
    template<class faceType1, class faceType2>
    inline label sharedVertex(const faceType1& f1, const faceType2& f2);

    //- find a vertex label common to all faces in the list
    template<class faceListType>
    inline label sharedVertex(const faceListType& fcs);

    //- check if selected elements are in one singly-connected chain
    template<class boolListType>
    inline bool areElementsInChain(const boolListType& sel);

    //- sort the edge chain such that the end of one edge is the beginning
    //- of the next one (valid for singly-connected chains, only)
    inline labelList sortEdgeChain(const DynList<edge>& bEdges);

    //- creates closed edge chains from the open chain
    inline void zipOpenChain(DynList<edge>& bEdges);

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "helperFunctionsTopologyManipulationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
