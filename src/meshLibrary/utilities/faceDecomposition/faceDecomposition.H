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
    faceDecomposition

Description
    Decomposition of a face into triangles

SourceFiles
    faceDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef faceDecomposition_H
#define faceDecomposition_H

#include "faceList.H"
#include "pointField.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faceDecomposition Declaration
\*---------------------------------------------------------------------------*/

class faceDecomposition
{
    // private data
        const face& f_;
    
        const pointField& points_;

    // private member functions
        //- find concave vertex and return its position
        //- in the face
        label concaveVertex() const;
    
        //- decomposes the face into triangle starting from
        //- the given vertex
        faceList decomposeFaceIntoTriangles(const label cv) const;

        //- disallows bitwise construct
        void operator=(const faceDecomposition&);

        //- copy constructor
        faceDecomposition(const faceDecomposition&);

public:

    // Constructors

        //- construct components
        faceDecomposition
        (
            const face&,
            const pointField&
        );

        //- Destructor
        ~faceDecomposition();

    // Member functions
        //- check if the face is convex
        bool isFaceConvex() const;

        //- check if the face is planar
        bool isFacePlanar() const;
        bool isFacePlanar(const scalar tol) const;

        //- decompose face into triangles
        faceList decomposeFaceIntoTriangles() const;
        
        //- decompose face into the minimal number
        //- of convex faces
        faceList decomposeFace() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
