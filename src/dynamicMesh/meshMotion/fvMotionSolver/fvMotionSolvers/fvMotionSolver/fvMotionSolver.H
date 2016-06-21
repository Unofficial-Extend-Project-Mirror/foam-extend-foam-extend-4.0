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
    Foam::fvMotionSolver

Description
    Virtual base class for finite volume mesh motion solvers.

    The boundary motion is set as a boundary condition on the motion velocity
    variable motionU.

SourceFiles
    fvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fvMotionSolver_H
#define fvMotionSolver_H

#include "motionSolver.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class fvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class fvMotionSolver
:
    public motionSolver
{

protected:

    // Protected data

        //- The fvMesh to be moved
        const fvMesh& fvMesh_;


    // Protected member functions

        //- Create the corresponding patch types for cellMotion from those
        //  of the given pointMotion
        template<class Type>
        wordList cellMotionBoundaryTypes
        (
            const typename GeometricField<Type, pointPatchField, pointMesh>::
            GeometricBoundaryField& pmUbf
        ) const;


public:

    //- Runtime type information
    TypeName("fvMotionSolver");


    // Constructors

        //- Construct from polyMesh
        fvMotionSolver(const polyMesh& mesh);


    // Destructor

        virtual ~fvMotionSolver();


    // Member Functions

        //- Return reference to the fvMesh to be moved
        const fvMesh& mesh() const
        {
            return fvMesh_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const = 0;

        //- Solve for motion
        virtual void solve() = 0;

        //- Update the pointMesh corresponding to the new points
        void movePoints(const pointField&);

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvMotionSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
