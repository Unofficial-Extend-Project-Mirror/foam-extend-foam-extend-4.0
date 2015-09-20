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

#include "thoboisMesh.H"
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
#include "tetMotionSolver.H"

#include "fixedValueTetPolyPatchFields.H"
#include "mixedTetPolyPatchFields.H"
#include "slipTetPolyPatchFields.H"
#include "zeroGradientTetPolyPatchFields.H"



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::thoboisMesh::makeLayersLive()
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



void Foam::thoboisMesh::valveDetach()
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

void Foam::thoboisMesh::valveAttach()
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

void Foam::thoboisMesh::prepareValveDetach()
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


bool Foam::thoboisMesh::update()
{
    Info << "bool Foam::layerSmooth::update()" << endl;
    tetMotionSolver& mSolver =
        refCast<tetMotionSolver>(msPtr_());

    tetPointVectorField& motionU = mSolver.motionU();

    // Find piston mesh modifier

    const label pistonLayerID =
        topoChanger_.findModifierID("pistonLayer");

    {
        valveDetach();

        topoChanger_[pistonLayerID].disable();

        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

        if (topoChangeMap->morphing())
        {
            mSolver.updateMesh(topoChangeMap());
            Info << "mSolver.updateMesh(topoChangeMap())" << endl;
        }

    }

    scalar deltaZ = engTime().pistonDisplacement().value();

    // deltaZ set to zero, FIXED PISTON POSITION
    deltaZ = 0.0;

    Info<< "deltaZ = " << deltaZ << " Piston at:" << pistonPosition()
        << endl;
    virtualPistonPosition() += deltaZ;

    Info << "pistonLayerID: " << pistonLayerID << endl;

    bool realDeformation = deformation();

    if
    (
        virtualPistonPosition() + deltaZ
      > deckHeight() - engTime().clearance().value() - SMALL
    )
    {
        realDeformation = true;
    }

    if (realDeformation)
    {
        // Disable layer addition
        Info << "**Mesh deformation mode" << endl;
        topoChanger_[pistonLayerID].disable();
    }
    else
    {
        // Enable layer addition
        Info << "**Piston layering mode" << endl;
        topoChanger_[pistonLayerID].enable();
    }

    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    bool meshChanged = topoChangeMap->morphing();

    pointField newPoints = allPoints();

    if (meshChanged)
    {
        mSolver.updateMesh(topoChangeMap());

        if (topoChangeMap().hasMotionPoints())
        {
            movePoints(topoChangeMap().preMotionPoints());
            newPoints = topoChangeMap().preMotionPoints();
        }
        setV0();
        resetMotion();
    }

    if (!deformation())
    {
//#       include "movePistonPointsLayering.H"
        movePoints(newPoints);

#       include "setValveMotionBoundaryConditionThobois.H"

        // Set piston velocity
        if (piston().patchID().active())
        {

            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[piston().patchID().index()]
                );

            pp.refValue() = vector::zero;
        }

        motionU.correctBoundaryConditions();


        DynamicList<label> constrainedPoints(mSolver.curPoints()().size()/100);
        DynamicList<vector> constrainedVelocity
        (
            mSolver.curPoints()().size()/100
        );

#       include "setThoboisMeshConstraintsNoDeformation.H"


        // Set piston velocity
        if (piston().patchID().active())
        {
            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[piston().patchID().index()]
                );

            pp.refValue() = vector::zero;
        }

        motionU.correctBoundaryConditions();

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

//         mSolver.solve
//         (
//             labelList(constrainedPoints.shrink()),
//             constrainedVelocity.shrink()
//         );

        mSolver.solve();

        //set to zero the motion U along the x and y directions

        newPoints = mSolver.curPoints();
        movePoints(newPoints);
        setVirtualPositions();
        mSolver.clearConstraints();
    }
    else
    {
#       include "setValveMotionBoundaryConditionThobois.H"
#       include "setPistonMotionBoundaryConditionThobois.H"
        DynamicList<label> constrainedPoints(mSolver.curPoints()().size()/100);
        DynamicList<vector> constrainedVelocity
        (
            mSolver.curPoints()().size()/100
        );

#       include "setThoboisMeshConstraints.H"

        // Set piston velocity
        if (piston().patchID().active())
        {

            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[piston().patchID().index()]
                );

            pp.refValue() = vector::zero;
        }

//        motionU.correctBoundaryConditions();

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

//        mSolver.solve
//        (
//            labelList(constrainedPoints.shrink()),
//            constrainedVelocity.shrink()
//        );

        mSolver.solve();

        //set to zero the motion U along the x and y directions

        newPoints = mSolver.curPoints();
        movePoints(newPoints);
        setVirtualPositions();
        mSolver.clearConstraints();
    }

    {
        pointField oldPointsNew = oldPoints();
        pointField newPointsNew = points();

        prepareValveDetach();
        topoChanger_[pistonLayerID].disable();
        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

//        Info << motionU << endl;

        if (topoChangeMap->morphing())
        {
            mSolver.updateMesh(topoChangeMap());

           {
               // correct the motion after attaching the sliding interface

               pointField mappedOldPointsNew(allPoints().size());

               mappedOldPointsNew.map(oldPointsNew, topoChangeMap->pointMap());
               pointField newPoints = allPoints();

               movePoints(mappedOldPointsNew);
               resetMotion();
               setV0();
               movePoints(newPoints);
           }
        }
    }

    Info << "CHANGED LAST" << endl;

/*
    if(moving() && meshChanged)
    {
        Info << "MOOOOOOOOOOVING!!!!!!!" << endl;
        Info << "min V0() post motion = " << min(V0()) << endl;
    }

    Info << "min V() post-motion = " << min(V()) << endl;
    Info << "max phi() post-motion = " << max(phi()) << endl;
    Info << "min phi() post-motion = " << min(phi()) << endl;
*/

    Info<< "Total cylinder volume at CA " << engTime().timeName() << " = "
        << sum(V()) << endl;

    return meshChanged;
}
