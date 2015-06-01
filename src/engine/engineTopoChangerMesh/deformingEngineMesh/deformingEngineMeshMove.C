/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

\*---------------------------------------------------------------------------*/

#include "deformingEngineMesh.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"

#include "tetPolyMesh.H"
#include "tetPointFields.H"
#include "elementFields.H"
#include "fixedValueTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"
#include "symmetryTetPolyPatch.H"

#include "tetFem.H"

#include "symmetryFvPatch.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"
#include "zeroGradientTetPolyPatchFields.H"
#include "tetMotionSolver.H"

#include "fixedValueTetPolyPatchFields.H"
#include "mixedTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"
#include "zeroGradientTetPolyPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::deformingEngineMesh::update()
{
    tetMotionSolver& mSolver =
        refCast<tetMotionSolver>(msPtr_());

    scalar deltaZ = engTime().pistonDisplacement().value();

    // deltaZ set to zero, FIXED PISTON POSITION
    deltaZ = 0.0;

    virtualPistonPosition() += deltaZ;

    pointField newPoints = points();

    {
#       include "setValveMotionBoundaryConditionDeformingEngineMesh.H"
#       include "setPistonMotionBoundaryConditionDeformingEngineMesh.H"
        Info << "piston motion" << endl;

        DynamicList<label> constrainedPoints(mSolver.curPoints()().size()/100);
        DynamicList<vector> constrainedVelocity
        (
            mSolver.curPoints()().size()/100
        );

#       include "setDeformingEngineMeshConstraints.H"


        labelList constrainedPointsList(constrainedPoints.shrink());
        vectorField constrainedVelocityField(constrainedVelocity.shrink());

        forAll (constrainedPointsList, i)
        {
            mSolver.setConstraint
            (
                constrainedPointsList[i],
                constrainedVelocityField[i]
            );
        }

        mSolver.solve();

        newPoints = mSolver.curPoints();
        movePoints(newPoints);
        setVirtualPositions();
        mSolver.clearConstraints();

    }

    pointField oldPointsNew = oldAllPoints();
    newPoints = points();
    movePoints(newPoints);
    Info << "max mesh phi before = " << max((phi())) << endl;
    Info << "min mesh phi before = " << min((phi())) << endl;


    Info << "max mesh phi after = " << max((phi())) << endl;
    Info << "min mesh phi after = " << min((phi())) << endl;

    return false;
}
