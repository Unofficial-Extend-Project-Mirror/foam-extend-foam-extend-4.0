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

Author
    Frank Bos, TU Delft.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "RBFMotionFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "RBFMotionSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(RBFMotionFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        RBFMotionFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBFMotionFunctionObject::RBFMotionFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    rotationAmplitude_(readScalar(dict.lookup("rotationAmplitude"))),
    rotationFrequency_(readScalar(dict.lookup("rotationFrequency"))),
    translationAmplitude_(dict.lookup("translationAmplitude")),
    translationFrequency_(dict.lookup("translationFrequency")),
    initialRotationOrigin_(dict.lookup("initialRotationOrigin")),
    statPoints_()

{
    Info << "Creating RBFMotion function object" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RBFMotionFunctionObject::start()
{
    const polyMesh& mesh =
        time_.lookupObject<polyMesh>(regionName_);

    // Grab RBF motion solver
    RBFMotionSolver& ms =
        const_cast<RBFMotionSolver&>
        (
            mesh.lookupObject<RBFMotionSolver>("dynamicMeshDict")
        );

    statPoints_ = ms.movingPoints();

    return true;
}


bool Foam::RBFMotionFunctionObject::execute()
{
    const polyMesh& mesh =
        time_.lookupObject<polyMesh>(regionName_);

    // Grab RBF motion solver
    RBFMotionSolver& ms =
        const_cast<RBFMotionSolver&>
        (
            mesh.lookupObject<RBFMotionSolver>("dynamicMeshDict")
        );

#   include "kinematicModel.H"

    ms.setMotion(motion);
    movePoints(ms.newPoints());

    return true;
}


bool Foam::RBFMotionFunctionObject::read(const dictionary& dict)
{
    return false;
}


// ************************************************************************* //
