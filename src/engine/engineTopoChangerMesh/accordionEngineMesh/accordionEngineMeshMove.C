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

\*---------------------------------------------------------------------------*/

#include "accordionEngineMesh.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "attachDetach.H"
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
#include "tetDecompositionMotionSolver.H"

#include "fixedValueTetPolyPatchFields.H"
#include "mixedTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"
#include "zeroGradientTetPolyPatchFields.H"



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::accordionEngineMesh::makeLayersLive()
{
    // Enable layering
    forAll (topoChanger_, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanger_[modI]))
        {
            topoChanger_[modI].enable();
        }
        else if (isA<slidingInterface>(topoChanger_[modI]))
        {
            topoChanger_[modI].disable();
        }
        else if (isA<attachDetach>(topoChanger_[modI]))
        {
            topoChanger_[modI].enable();
        }
        else
        {
            FatalErrorIn("void Foam::engineTopoFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanger_[modI].type()
                << abort(FatalError);
        }
    }
}

void Foam::accordionEngineMesh::makeDetachLive()
{
    // Enable sliding interface
    forAll (topoChanger_, modI)
    {
        if (isA<attachDetach>(topoChanger_[modI]))
        {
            topoChanger_[modI].enable();
        }
        else
        {
            FatalErrorIn("void Foam::accordionEngineMesh::makeDetachLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanger_[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::accordionEngineMesh::valveDetach()
{
    // Enable sliding interface
    forAll (topoChanger_, modI)
    {
        if (isA<attachDetach>(topoChanger_[modI]))
        {
            const attachDetach& ad =
                refCast<const attachDetach>(topoChanger_[modI]);

            const word masterName = ad.masterPatchID().name();

            // Find the valve with that name
            label valveIndex = -1;

            forAll (valves_, valveI)
            {
                if
                (
                    valves_[valveI].detachInCylinderPatchID().name()
                 == masterName
                )
                {
                    valveIndex = valveI;
                    break;
                }
            }

            if (valveIndex < 0)
            {
                FatalErrorIn
                (
                    "void Foam::engineTopoFvMesh::prepareValveDetach()"
                )   << "Cannot match patch for attach/detach " << modI
                    << abort(FatalError);
            }

            if (debug)
            {
                Info<< " valveI: " << valveIndex << " attached: "
                    << ad.attached()
                    << " valve open: " << valves_[valveIndex].isOpen()
                    << endl;
            }

            ad.setDetach();
        }
    }
}

void Foam::accordionEngineMesh::valveAttach()
{
    // Enable sliding interface
    forAll (topoChanger_, modI)
    {
        if (isA<attachDetach>(topoChanger_[modI]))
        {
            const attachDetach& ad =
                refCast<const attachDetach>(topoChanger_[modI]);

            const word masterName = ad.masterPatchID().name();

            // Find the valve with that name
            label valveIndex = -1;

            forAll (valves_, valveI)
            {
                if
                (
                    valves_[valveI].detachInCylinderPatchID().name()
                 == masterName
                )
                {
                    valveIndex = valveI;
                    break;
                }
            }

            if (valveIndex < 0)
            {
                FatalErrorIn
                (
                    "void Foam::engineTopoFvMesh::prepareValveDetach()"
                )   << "Cannot match patch for attach/detach " << modI
                    << abort(FatalError);
            }

            if (debug)
            {
                Info<< " valveI: " << valveIndex << " attached: "
                    << ad.attached()
                    << " valve open: " << valves_[valveIndex].isOpen()
                    << endl;
            }

            ad.setAttach();
        }
    }
}


void Foam::accordionEngineMesh::prepareValveDetach()
{
    // Enable sliding interface
    forAll (topoChanger_, modI)
    {
        if (isA<attachDetach>(topoChanger_[modI]))
        {
            const attachDetach& ad =
                refCast<const attachDetach>(topoChanger_[modI]);

            const word masterName = ad.masterPatchID().name();

            // Find the valve with that name
            label valveIndex = -1;

            forAll (valves_, valveI)
            {
                if
                (
                    valves_[valveI].detachInCylinderPatchID().name()
                 == masterName
                )
                {
                    valveIndex = valveI;
                    break;
                }
            }

            if (valveIndex < 0)
            {
                FatalErrorIn
                (
                    "void Foam::engineTopoFvMesh::prepareValveDetach()"
                )   << "Cannot match patch for attach/detach " << modI
                    << abort(FatalError);
            }

            if (debug)
            {
                Info<< " valveI: " << valveIndex << " attached: "
                    << ad.attached()
                    << " valve open: " << valves_[valveIndex].isOpen()
                    << endl;
            }

            if (valves_[valveIndex].isOpen())
            {
                ad.setAttach();
            }
            else
            {
                ad.setDetach();
            }
        }
    }
}


bool Foam::accordionEngineMesh::update()
{
    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

    tetPointVectorField& motionU = mSolver.motionU();

    Info << "motioU.size() = " << motionU.internalField().size() << endl;

    scalar deltaZ = engTime().pistonDisplacement().value();

    // deltaZ set to zero, FIXED PISTON POSITION
    deltaZ = 0.0;

    virtualPistonPosition() += deltaZ;

    makeDetachLive();

//    valveDetach();

    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    if (topoChangeMap->morphing())
    {
        mSolver.updateMesh(topoChangeMap());
        Info << "mSolver.updateMesh(topoChangeMap())" << endl;

        if (topoChangeMap->hasMotionPoints())
        {
            movePoints(topoChangeMap->preMotionPoints());
            resetMotion();
            setV0();
        }
    }

    pointField newPoints = points();

    {
#       include "setValveMotionBoundaryConditionAccordionEngineMesh.H"
#       include "setPistonMotionBoundaryConditionAccordionEngineMesh.H"


        DynamicList<label> constrainedPoints(mSolver.curPoints()().size()/100);
        DynamicList<vector> constrainedVelocity
        (
            mSolver.curPoints()().size()/100
        );

#       include "setAccordionEngineMeshConstraints.H"


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

        //set to zero the motion U along the x and y directions

        newPoints = mSolver.curPoints();
        movePoints(newPoints);
        setVirtualPositions();
        mSolver.clearConstraints();
    }

//    pointField oldPointsNew = oldPoints();
    pointField oldPointsNew = oldAllPoints();
    newPoints = points();

    Info << "max mesh phi before = " << max((phi())) << endl;
    Info << "min mesh phi before = " << min((phi())) << endl;

//    makeDetachLive();
    prepareValveDetach();

    autoPtr<mapPolyMesh> topoChangeMap1 = topoChanger_.changeMesh();

    if (topoChangeMap1->morphing())
    {
        mSolver.updateMesh(topoChangeMap1());

        if (topoChangeMap1->hasMotionPoints())
        {
//              movePoints(topoChangeMap1->preMotionPoints());
//              resetMotion();
//              setV0();

            // correct the motion after attaching the sliding interface
            pointField mappedOldPointsNew(allPoints().size());

            mappedOldPointsNew.map(oldPointsNew, topoChangeMap1->pointMap());
            pointField newPoints = allPoints();

            movePoints(mappedOldPointsNew);
            resetMotion();
            setV0();
            movePoints(newPoints);

        }
    }

    Info << "max mesh phi after = " << max((phi())) << endl;
    Info << "min mesh phi after = " << min((phi())) << endl;

/*         if(correctPointsMotion_)
        {
            // correct the motion after attaching the sliding interface
            pointField mappedOldPointsNew(allPoints().size());

            mappedOldPointsNew.map(oldPointsNew, topoChangeMap3->pointMap());
            pointField newPoints = allPoints();

            movePoints(mappedOldPointsNew);
            resetMotion();
            setV0();
            movePoints(newPoints);
        }

//    Info << motionU << endl;

    if (meshChanged1)
    {
            Info << "meshChanged1" << endl;

        mSolver.updateMesh(topoChangeMap1());

       {
           pointField mappedOldPointsNew(allPoints().size());

           mappedOldPointsNew.map(oldPointsNew, topoChangeMap1->pointMap());
           pointField newPoints = allPoints();

           movePoints(mappedOldPointsNew);
           resetMotion();
           setV0();
           movePoints(newPoints);
       }
    }
*/
    Info<< "Total cylinder volume at CA " << engTime().timeName() << " = "
        << sum(V()) << endl;

    return false;
}
