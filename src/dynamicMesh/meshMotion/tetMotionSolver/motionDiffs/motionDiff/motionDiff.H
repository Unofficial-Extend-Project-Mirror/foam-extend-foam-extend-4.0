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
