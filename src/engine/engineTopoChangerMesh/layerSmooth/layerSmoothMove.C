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

#include "layerSmooth.H"
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

void Foam::layerSmooth::makeLayersLive()
{
    // Enable layering
    forAll (topoChanger_, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanger_[modI]))
        {
            topoChanger_[modI].enable();
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


void Foam::layerSmooth::valveDetach()
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


void Foam::layerSmooth::valveAttach()
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


void Foam::layerSmooth::prepareValveDetach()
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


bool Foam::layerSmooth::update()
{

    Info << "bool Foam::layerSmooth::update()" << endl;

    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

    tetPointVectorField& motionU = mSolver.motionU();

//    motionU.internalField() = (vector::zero);

//    Info << motionU << endl;


    // Find piston mesh modifier
    const label pistonLayerID = topoChanger_.findModifierID("pistonLayer");

    {
        valveDetach();

        Info << "valveDetach()" << endl;

        topoChanger_[pistonLayerID].disable();

        Info << "disable()" << endl;

        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

        Info << "changeMesh()" << endl;

        if (topoChangeMap->morphing())
        {
            mSolver.updateMesh(topoChangeMap());
            Info << "mSolver.updateMesh(topoChangeMap())" << endl;
        }

    }

    scalar deltaZ = engTime().pistonDisplacement().value();
    Info<< "deltaZ = " << deltaZ << " Piston at:" << pistonPosition()
        << endl;
    virtualPistonPosition() += deltaZ;

    Info << "pistonLayerID: " << pistonLayerID << endl;

    const layerAdditionRemoval& pistonLayers =
        dynamicCast<const layerAdditionRemoval>
        (
            topoChanger_[pistonLayerID]
        );

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

    scalar minLayerThickness = pistonLayers.minLayerThickness();

    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    pointField newPoints = points();

    if (topoChangeMap->morphing())
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

    if(!deformation())
    {
#       include "movePistonPointsLayeringLayerSmooth.H"
        Info << "movePoints" << endl;
        movePoints(newPoints);
        Info << "setBoundaryMotion" << endl;

//#       include "setBoundaryMotion.H"

#       include "setValveMotionBoundaryCondition.H"

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


//#       include "setPistonMotionBoundaryCondition.H"
        Info << "mSolver" << endl;
        mSolver.solve();
        Info << "newPoints" << endl;
        newPoints = mSolver.curPoints();
        Info << "movePoints 2" << endl;
        movePoints(newPoints);
    }
    else
    {
#       include "setValveMotionBoundaryCondition.H"
#       include "setPistonMotionBoundaryConditionLayerSmooth.H"
        mSolver.solve();
        newPoints = mSolver.curPoints();
        movePoints(newPoints);
    }

    {
        pointField oldPointsNew = oldPoints();
        pointField newPointsNew = points();

        prepareValveDetach();
        topoChanger_[pistonLayerID].disable();
        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

        Info << "changeMesh" << endl;


        if (topoChangeMap->morphing())
        {
            Info << "meshChanged" << endl;

            mSolver.updateMesh(topoChangeMap());
            Info << "updateMesh" << endl;

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

//    Info << motionU << endl;

#   ifdef CheckMesh
    checkMesh(true);
#   endif

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

    Info << "Total cylinder volume at CA " << engTime().timeName() << " = " <<
        sum(V()) << endl;

    return meshChanged;
}
/*

   {

       pointField oldPointsNew = oldPoints();
       pointField newPointsNew = points();

       // Attach the interface
       Info << "Coupling sliding interfaces" << endl;
       makeSlidersLive();
       valveAttach();
       prepareValveDetach();



       // Changing topology by hand
       autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();
       Info << "Sliding interfaces coupled: " << attached() << endl;

       if (topoChangeMap3->morphing())
       {
           mSolver.updateMesh(topoChangeMap3());

           if (topoChangeMap3->hasMotionPoints())
           {
//                movePoints(topoChangeMap3->preMotionPoints());
//                resetMotion();
//                setV0();
           }
                                       if(correctPointsMotion_)
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
       }
     }

*/
