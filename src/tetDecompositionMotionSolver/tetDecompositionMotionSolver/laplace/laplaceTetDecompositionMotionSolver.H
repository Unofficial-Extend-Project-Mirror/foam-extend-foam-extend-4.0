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
    laplaceTetDecompositionMotionSolver

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based Laplace motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

SourceFiles
    laplaceTetDecompositionMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef laplaceTetDecompositionMotionSolver_H
#define laplaceTetDecompositionMotionSolver_H

#include "tetDecompositionMotionSolver.H"
#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "lduMatrix.H"
#include "motionDiff.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionDiff;

/*---------------------------------------------------------------------------*\
                           Class laplaceTetDecompositionMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class laplaceTetDecompositionMotionSolver
:
    public tetDecompositionMotionSolver
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        laplaceTetDecompositionMotionSolver
        (
            const laplaceTetDecompositionMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const laplaceTetDecompositionMotionSolver&);


protected:

        //- Diffusion
        motionDiff* diffusionPtr_;

        //- Is diffusion field frozen
        bool frozen_;

        //- First motion
        bool firstMotion_;

        //- Return pointer to motion diffusivity
        const motionDiff& diffusion()
        {
            return *diffusionPtr_;
        }

        // Holds recent solver performance
        lduMatrix::solverPerformance solverPerf_;


public:

    //- Runtime type information
#ifdef CELL_DECOMP
    TypeName("laplaceCellDecomposition");
#elif defined(FACE_DECOMP)
    TypeName("laplaceFaceDecomposition");
#else
#error Undefined: cell or face decomposition
#endif


    // Constructors

        laplaceTetDecompositionMotionSolver
        (
            const polyMesh& mesh,
            Istream& msData
        );


    // Destructor

        virtual ~laplaceTetDecompositionMotionSolver();


    // Member Functions

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);

        //- Return recent solver performance
        const lduMatrix::solverPerformance& solverPerformance() const
        {
            return solverPerf_;
        }

        //- Return pointer to motion diffusivity
        const elementScalarField& motionGamma() const
        {
            return diffusionPtr_->motionGamma();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
