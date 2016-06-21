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
    subsetMotionSolverFvMesh

Description

SourceFiles
    subsetMotionSolverFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef subsetMotionSolverFvMesh_H
#define subsetMotionSolverFvMesh_H

#include "dynamicFvMesh.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class subsetMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class subsetMotionSolverFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Moving mesh dictionary
        dictionary movingMeshCoeffs_;

        //- Subset mesh
        fvMeshSubset subsetMesh_;

        //- Pointer to mesh motion
        autoPtr<motionSolver> motionPtr_;

        //- Motion under-relaxation coefficient
        scalar alpha_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        subsetMotionSolverFvMesh(const subsetMotionSolverFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const subsetMotionSolverFvMesh&);


public:

    //- Runtime type information
    TypeName("subsetMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        subsetMotionSolverFvMesh(const IOobject& io);


    // Destructor

        virtual ~subsetMotionSolverFvMesh();


    // Member Functions

        //- Return reference to subset mesh
        const fvMeshSubset& subsetMesh() const;

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
