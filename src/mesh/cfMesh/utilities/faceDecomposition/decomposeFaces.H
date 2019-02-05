/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    decomposeFaces

Description
    Decompose selected mesh faces into triangles using midnode subdivision

SourceFiles
    decomposeFaces.C

\*---------------------------------------------------------------------------*/

#ifndef decomposeFaces_H
#define decomposeFaces_H

#include "polyMeshGenModifier.H"
#include "DynList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class decomposeFaces Declaration
\*---------------------------------------------------------------------------*/

class decomposeFaces
{
    // private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- number of points
        VRWGraph newFacesForFace_;

        //- is decomposition performed
        bool done_;

        //- disallows bitwise construct
        void operator=(const decomposeFaces&);

        //- copy constructor
        decomposeFaces(const decomposeFaces&);

public:

    // Constructors

        //- construct components
        decomposeFaces(polyMeshGen& mesh);

        //- Destructor
        ~decomposeFaces();

    // Member functions
        //- decompose selected faces into triangles using midnode subdivision
        void decomposeMeshFaces(const boolList& decomposeFace);

        //- decompose internal faces containing concave nodes
        void decomposeConcaveInternalFaces(const boolList& concaveVertex);

        //- reference to a graph containg label of new faces generated
        //- from each face
        const VRWGraph& newFacesForFace() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
