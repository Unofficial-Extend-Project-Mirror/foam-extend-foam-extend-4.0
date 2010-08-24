/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2010 Tommaso Lucchini
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
    verticalValvesGambit

\*---------------------------------------------------------------------------*/

#include "pistonSliding.H"
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

#include "zeroGradientFvPatchFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pistonSliding::makeLayersLive()
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

void Foam::pistonSliding::makeSlidersLive()
{
    // Enable sliding interface
    forAll (topoChanger_, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanger_[modI]))
        {
            topoChanger_[modI].disable();
        }
        else if (isA<slidingInterface>(topoChanger_[modI]))
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

bool Foam::pistonSliding::attached() const
{
    const polyTopoChanger& morphs = topoChanger_;

    bool result = false;

    forAll (morphs, modI)
    {
        if (typeid(morphs[modI]) == typeid(slidingInterface))
        {
            result =
                result
             || refCast<const slidingInterface>(morphs[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll (morphs, modI)
    {
        if (typeid(morphs[modI]) == typeid(slidingInterface))
        {
            if
            (
                result
             != refCast<const slidingInterface>(morphs[modI]).attached()
            )
            {
                FatalErrorIn("bool movingSquaresTM::attached() const")
                    << "Slider " << modI << " named " << morphs[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    return result;
}


void Foam::pistonSliding::valveDetach()
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

void Foam::pistonSliding::valveAttach()
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

void Foam::pistonSliding::prepareValveDetach()
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


bool Foam::pistonSliding::update()
{

    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

     tetPointVectorField& motionU = mSolver.motionU();

    // Detaching the interfacethobois2DSlidingDeform
    if (attached())
    {
        Info << "Decoupling sliding interfaces" << endl;

        makeSlidersLive();
        valveDetach();
        autoPtr<mapPolyMesh> topoChangeMap1 = topoChanger_.changeMesh();

        Info << "sliding interfaces successfully decoupled!!!" << endl;
        if (topoChangeMap1->morphing())
        {
            mSolver.updateMesh(topoChangeMap1());
        }
    }
    else
    {
        Info << "Sliding interfaces decoupled" << endl;
        valveDetach();
    }

    Info << "Executing layer action" << endl;

    // Layering mode
    makeLayersLive();

    const label pistonLayerID =
        topoChanger_.findModifierID("pistonLayer");

    scalar deltaZ = engTime().pistonDisplacement().value();

    if (deformation())
    {
        Info << "Mesh deformation mode" << endl;
        topoChanger_[pistonLayerID].disable();
    }
    else
    {
        // Activate piston layer
        Info << "Piston layering mode" << endl;
        topoChanger_[pistonLayerID].enable();
    }

    {
        // Activate bottom layer
        Info << "Valve layering mode" << endl;
        forAll(valves_, valveI)
        {
            if(valves_[valveI].bottomPatchID().active())
            {
                // Find piston mesh modifier
                const label valveLayerID2 =
                    topoChanger_.findModifierID
                    (
                        "valveBottomLayer"
                      + Foam::name(valveI + 1)
                    );
                topoChanger_[valveLayerID2].enable();
            }
        }
    }

    forAll(valves_, valveI)
    {
        if(!realDeformation())
        {
            scalar valveDisplacement =
                valves_[valveI].curVelocity()*
                valves_[valveI].cs().axis().z()*engTime().deltaT().value();

            Info<< "valveDisplacement = " << valveDisplacement << nl
                << "valve velocity =" << valves_[valveI].curVelocity() << endl;

            valveTopPosition_[valveI] += valveDisplacement;
            valveBottomPosition_[valveI] += valveDisplacement;
            valvePistonPosition_[valveI] += deltaZ;
        }
    }

    forAll(valves_,valveI)
    {
        if(valves_[valveI].curLift() > valves_[valveI].deformationLift())
        {
            if(valves_[valveI].poppetPatchID().active())
            {
                // Find valve layer mesh modifier
                const label valveLayerID =
                    topoChanger_.findModifierID
                    (
                        "valvePoppetLayer"
                      + Foam::name(valveI + 1)
                    );

                if (valveLayerID < 0)
                {
                    FatalErrorIn("void engineFvMesh::moveAndMorph()")
                        << "valve modifier not found."
                        << abort(FatalError);
                }

                topoChanger_[valveLayerID].enable();
            }
        }
        else
        {
            if(valves_[valveI].poppetPatchID().active())
            {
                // Find valve layer mesh modifier
                const label valveLayerID =
                    topoChanger_.findModifierID
                    (
                        "valvePoppetLayer"
                      + Foam::name(valveI + 1)
                    );

                if (valveLayerID < 0)
                {
                    FatalErrorIn("void engineFvMesh::moveAndMorph()")
                        << "valve modifier not found."
                        << abort(FatalError);
                }

                topoChanger_[valveLayerID].disable();
            }
        }
    }

    pointField newPoints = points();

//    pointField newPoints = allPoints();

    autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

    // Changing topology by hand
    if (topoChangeMap2->morphing())
    {

        mSolver.updateMesh(topoChangeMap2());

        if (topoChangeMap2->hasMotionPoints())
        {
            movePoints(topoChangeMap2->preMotionPoints());
            newPoints = topoChangeMap2->preMotionPoints();
        }
        setV0();
        resetMotion();

    }

    // Work array for new points position.

// NEW !

    // Reset the position of layered interfaces

    {

#        include "setValveMotionBoundaryConditionPistonSliding.H"

        DynamicList<label> constrainedPoints(mSolver.curPoints()().size()/100);
        DynamicList<vector> constrainedVelocity
        (
            mSolver.curPoints()().size()/100
        );

#       include "setPistonSlidingConstraint.H"

        forAll (constrainedPoints, i)
        {
               mSolver.setConstraint
               (
                   constrainedPoints[i],
                   constrainedVelocity[i]
               );
        }

//        mSolver.solve();

        newPoints = mSolver.curPoints();
        movePoints(newPoints);

        setVirtualPositions();
        mSolver.clearConstraints();

//      layering

#       include "moveValvePointsPistonSliding.H"
//#       include "movePistonPointsPistonSliding.H"
        movePoints(newPoints);
//        newPoints = points();
        newPoints = allPoints();
        setVirtualPositions();

    }

    forAll(valves(), valveI)
    {
        Info << "Valve Top Position for valve n. " << valveI + 1 << " = " <<
        valveTopPosition_[valveI] << endl;
        Info << "Valve Bottom Position for valve n. " << valveI + 1 << " = " <<
        valveBottomPosition_[valveI] << endl;
    }

    Info << "virtualPistonPosition = " << virtualPistonPosition()
    << ", deckHeight = " << deckHeight()
    << ", pistonPosition = " << pistonPosition()
    << endl;

    pistonPosition() += deltaZ;

    scalarField Vold = V();

    volScalarField vMesh
    (
        IOobject
        (
            "vMesh",
            time().timeName(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    );

    vMesh.internalField() = V();
    vMesh.oldTime().internalField() = V0();


    // Coupling the interface Again

    {

        pointField newPointsNew = allPoints();

        // Attach the interface
        Info << "Coupling sliding interfaces" << endl;
        makeSlidersLive();
        valveAttach();
        prepareValveDetach();

        pointField oldPointsNew = oldAllPoints();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        Info << "Sliding interfaces coupled: " << attached() << endl;
        if (topoChangeMap3->morphing())
        {
            Info << "mesh changed 3" << endl;

            mSolver.updateMesh(topoChangeMap3());


            if (topoChangeMap3->hasMotionPoints())
            {
                movePoints(topoChangeMap3->preMotionPoints());
                resetMotion();
                setV0();
            }

            {
                // correct the motion after attaching the sliding interface
                Info << "oldPointsNew.size = " << oldPointsNew.size() << nl
                    << "allPointsNew.size = " << allPoints().size() << nl
                    << "pointsNew.size = " << points().size() << endl;

                pointField mappedOldPointsNew(allPoints().size());

                pointField newPoints = allPoints();
//                pointField newPoints = points();

                mappedOldPointsNew.map
                (
                    oldPointsNew,
                    topoChangeMap3->pointMap()
                );
                movePoints(mappedOldPointsNew);

                resetMotion();
                setV0();

                movePoints(mappedOldPointsNew);

                Info<< "mappedOldPointsNew.size() = "
                    << mappedOldPointsNew.size() << nl
                    << "newPoints.size() = " << newPoints.size() << nl
                    << "topoChangeMap3->pointMap().size() = "
                    <<  topoChangeMap3->pointMap().size() << endl;

                movePoints(newPoints);
            }
        }
    }

    if(moving())
    {
        Info << "Mesh moving, OK" << endl;
    }
    else
    {
        Info << "Mesh NOT moving, WARNINGGGG" << endl;
    }

    return true;
}


