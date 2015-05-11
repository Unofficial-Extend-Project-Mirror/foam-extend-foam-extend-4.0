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
    Foam::refineImmersedBoundaryMesh

Description
    Refine a mesh with immersed boundary

SourceFiles
    refineImmersedBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef refineImmersedBoundaryMesh_H
#define refineImmersedBoundaryMesh_H

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class refineImmersedBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class refineImmersedBoundaryMesh
{
public:

    // Public enumerations

        //- Cell collection method
        enum ibCellCollection
        {
            UNDEFINED_COLLECTION = 0,
            IB_CELLS,
            IB_CELL_CELLS,
            IB_CELL_CELL_FACES
        };

        //- Projection method names
        static const NamedEnum<ibCellCollection, 4> ibCellCollectionNames_;


private:

    // Private data

        //- Reference to mesh
        fvMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        refineImmersedBoundaryMesh(const refineImmersedBoundaryMesh&);

        //- Disallow default bitwise assignment
        void operator=(const refineImmersedBoundaryMesh&);


        // Mesh refinement

            //- Add ibCells for refinement
            void addIbCells(labelHashSet& refCellSet) const;

            //- Add ibCellCells for refinement
            void addIbCellCells(labelHashSet& refCellSet) const;

            //- Add ibCellCellFaces for refinement
            void addIbCellCellFaces(labelHashSet& refCellSet) const;


            //- From refineMesh application
            label axis(const vector& normal) const;

            //- From refineMesh application
            label twoDNess() const;


public:

    // Constructors

        //- Construct from mesh
        refineImmersedBoundaryMesh(fvMesh& mesh);


    //- Destructor
    ~refineImmersedBoundaryMesh();


    // Member Functions

        //- Return list of cells for refinement based on specifie collection
        //  type (see enumeration)
        labelList refinementCells
        (
            const ibCellCollection& collectionType
        ) const;

        // Refine mesh given a list of cells for refinement
        void refineMesh(const labelList& refCells) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
