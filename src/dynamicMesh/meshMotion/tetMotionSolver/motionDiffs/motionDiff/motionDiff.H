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
    motionDiff

Description
    Virtual base class for motion diffusion.

SourceFiles
    motionDiff.C

\*---------------------------------------------------------------------------*/

#ifndef motionDiff_H
#define motionDiff_H

#include "tetPolyMesh.H"
#include "elementFieldsFwd.H"
#include "tetMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class motionDiff Declaration
\*---------------------------------------------------------------------------*/

class motionDiff
{
private:

    // Private data

        //- Mesh reference
        const tetMotionSolver& mSolver_;


public:

    //- Runtime type information
    TypeName("motionDiff");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            motionDiff,
            dictionary,
            (
                const tetMotionSolver& mSolver
            ),
            (mSolver)
        );


    // Selectors

        //- Select null constructed
        static autoPtr<motionDiff> New
        (
            const tetMotionSolver& mSolver
        );


    // Constructors

        //- Construct from components
        motionDiff(const tetMotionSolver& mSolver)
        :
            mSolver_(mSolver)
        {}


    // Destructor

        virtual ~motionDiff();


    // Member Functions

        //- Return reference to mesh
        const tetMotionSolver& mSolver() const
        {
            return mSolver_;
        }

        //- Return reference to mesh
        const tetPolyMesh& tetMesh() const
        {
            return mSolver_.tetMesh();
        }

        //- Return diffusivity field
        virtual const elementScalarField& motionGamma() const = 0;


        //- Correct the motion diffusivity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
