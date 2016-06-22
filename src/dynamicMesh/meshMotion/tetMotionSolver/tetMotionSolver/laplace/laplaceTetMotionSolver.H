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
    laplaceTetMotionSolver

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based Laplace motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

SourceFiles
    laplaceTetMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef laplaceTetMotionSolver_H
#define laplaceTetMotionSolver_H

#include "tetMotionSolver.H"
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
                    Class laplaceTetMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class laplaceTetMotionSolver
:
    public tetMotionSolver
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        laplaceTetMotionSolver
        (
            const laplaceTetMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const laplaceTetMotionSolver&);


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
        lduSolverPerformance solverPerf_;


public:

    //- Runtime type information
    TypeName("laplace");


    // Constructors

        laplaceTetMotionSolver
        (
            const polyMesh& mesh,
            Istream& msData
        );


    // Destructor

        virtual ~laplaceTetMotionSolver();


    // Member Functions

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);

        //- Return recent solver performance
        const lduSolverPerformance& solverPerformance() const
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
