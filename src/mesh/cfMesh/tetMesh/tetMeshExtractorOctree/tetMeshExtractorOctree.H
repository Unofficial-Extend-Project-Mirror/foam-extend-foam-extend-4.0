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
    tetMeshExtractorOctree

Description
    A class which extracts tet mesh out of an octree structure

SourceFiles
    tetMeshExtractorOctree.C

\*---------------------------------------------------------------------------*/

#ifndef tetMeshExtractorOctree_H
#define tetMeshExtractorOctree_H

#include "polyMeshGenModifier.H"
#include "partTet.H"
#include "tetCreatorOctree.H"
#include "VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tetMeshExtractorOctree Declaration
\*---------------------------------------------------------------------------*/

class tetMeshExtractorOctree
{
    // Private data
    //- create tets
    tetCreatorOctree tetCreator_;

    //- reference to the mesh
    polyMeshGen& mesh_;

    // Private member functions
    //- copy tetPoints_ into mesh
    void createPoints();

    //- create mesh data
    void createPolyMesh();

    // Private copy constructor
    //- Disallow default bitwise copy construct
    tetMeshExtractorOctree ( const tetMeshExtractorOctree& );

    //- Disallow default bitwise assignment
    void operator= ( const tetMeshExtractorOctree& );

public:

    // Constructors

    //- Construct from octree and mesh data
    tetMeshExtractorOctree
    (
        const meshOctree& octree,
        const IOdictionary& meshDict,
        polyMeshGen& mesh
    );

    // Destructor

    ~tetMeshExtractorOctree();


    // Member Functions
    void createMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
