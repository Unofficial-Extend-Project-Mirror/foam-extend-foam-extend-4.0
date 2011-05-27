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
    RBFMotionSolver

Description
    Radial basis function motion solver

Author
    Frank Bos, TU Delft.  All rights reserved.

SourceFiles
    RBFMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef RBFMotionSolver_H
#define RBFMotionSolver_H

#include "motionSolver.H"
#include "polyMesh.H"
#include "RBFInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RBFMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class RBFMotionSolver
:
    public motionSolver
{
    // Private data

        //- Moving patch names
        wordList movingPatches_;

        //- Static patch names
        wordList staticPatches_;

        //- Coarsening ratio
        label coarseningRatio_;

        //- Include zero motion of static patches in RBF interpolation
        Switch includeStaticPatches_;

        //- Frozen interpolation
        Switch frozenInterpolation_;

        //- Moving point IDs
        labelList movingIDs_;

        //- Moving points on the boundary
        mutable vectorField movingPoints_;

        //- Static point IDs
        labelList staticIDs_;

        //- Control point IDs
        labelList controlIDs_;

        //- Control points on the boundary
        mutable vectorField controlPoints_;

        //- Internal point IDs
        labelList internalIDs_;

        //- Internal points
        vectorField internalPoints_;

        //- Prescribed motion on control points
        vectorField motion_;

        //- RBF interpolation
        RBFInterpolation interpolation_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        RBFMotionSolver(const RBFMotionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const RBFMotionSolver&);


        //- Make control point IDs.  Constructor helper
        void makeControlIDs();

        //- Set location of points
        void setMovingPoints() const;


public:

    //- Runtime type information
    TypeName("RBFMotionSolver");


    // Constructors

        //- Construct from polyMesh
        RBFMotionSolver
        (
            const polyMesh& mesh,
            Istream& msData
        );


    // Destructor

        virtual ~RBFMotionSolver();


    // Member Functions

        //- Return moving point indices in mesh points
        const labelList& movingIDs() const
        {
            return movingIDs_;
        }

        //- Return moving points
        const vectorField& movingPoints() const;

        //- Set motion of moving points
        void setMotion(const vectorField&);


        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
