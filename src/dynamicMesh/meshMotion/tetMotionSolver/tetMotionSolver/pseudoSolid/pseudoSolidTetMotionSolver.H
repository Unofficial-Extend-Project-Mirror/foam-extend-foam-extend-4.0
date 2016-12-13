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
    pseudoSolidTetMotionSolver

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based pseudo-solid motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

    Note: Diffusion field from laplaceTetDecompositionMotionSolver is used
          as Lame's constant \mu

SourceFiles
    pseudoSolidTetMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pseudoSolidTetMotionSolver_H
#define pseudoSolidTetMotionSolver_H

#include "laplaceTetMotionSolver.H"

#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class pseudoSolidTetMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class pseudoSolidTetMotionSolver
:
    public laplaceTetMotionSolver
{
    // Private data

        //- Poisson's ration
        scalar nu_;

        //- Number of correctors in solution procedure of pseudo-solid eqn
        label nCorrectors_;

        //- Convergence tolerance for solution of pseudo-solid eqn
        scalar convergenceTolerance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pseudoSolidTetMotionSolver
        (
            const pseudoSolidTetMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const pseudoSolidTetMotionSolver&);


public:

    //- Runtime type information
    TypeName("pseudoSolid");


    // Constructors

        pseudoSolidTetMotionSolver
        (
            const polyMesh& mesh,
            Istream& msData
        );


    // Destructor

        virtual ~pseudoSolidTetMotionSolver();


    // Member Functions

        //- Solve for motion
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
