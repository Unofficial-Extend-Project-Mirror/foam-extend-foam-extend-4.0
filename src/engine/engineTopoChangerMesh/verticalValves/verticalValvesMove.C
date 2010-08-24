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
    verticalValves

\*---------------------------------------------------------------------------*/

#include "verticalValves.H"
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

void Foam::verticalValves::makeLayersLive()
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

void Foam::verticalValves::makeSlidersLive()
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

bool Foam::verticalValves::attached() const
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

void Foam::verticalValves::valveDetach()
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

void Foam::verticalValves::valveAttach()
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

void Foam::verticalValves::prepareValveDetach()
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


bool Foam::verticalValves::update()
{

    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

    // Detaching the interface
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

    // Piston Layering

    makeLayersLive();
    scalar deltaZ = engTime().pistonDisplacement().value();

    // Find piston mesh modifier
    const label pistonLayerID =
        topoChanger_.findModifierID("pistonLayer");

    if (realDeformation())
    {
        Info << "Mesh deformation mode" << endl;
        topoChanger_[pistonLayerID].disable();
        forAll(valves_, valveI)
        {
            // Find piston mesh modifier
            const label valveLayerID1 =
                topoChanger_.findModifierID
                (
                    "valvePistonLayer"
                  + Foam::name(valveI + 1)
                );

            topoChanger_[valveLayerID1].disable();

            if(valves_[valveI].bottomPatchID().active())
            {
                // Find piston mesh modifier
                const label valveLayerID2 =
                    topoChanger_.findModifierID
                    (
                        "valveBottomLayer"
                      + Foam::name(valveI + 1)
                    );

                topoChanger_[valveLayerID2].disable();
            }
        }
    }
    else
    {
        // Activate piston layer
        Info << "Piston layering mode" << endl;
        topoChanger_[pistonLayerID].enable();
        forAll(valves_, valveI)
        {
            // Find piston mesh modifier
            const label valveLayerID1 =
                topoChanger_.findModifierID
                (
                    "valvePistonLayer"
                  + Foam::name(valveI + 1)
                );
            topoChanger_[valveLayerID1].enable();

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

    scalar minLayerThickness = piston().minLayer();
    virtualPistonPosition() += deltaZ;

    forAll(valves_,valveI)
    {
        if(!realDeformation())
        {
            scalar valveDisplacement =
                valves_[valveI].curVelocity()*
                valves_[valveI].cs().axis().z()*engTime().deltaT().value();

            valveTopPosition_[valveI] += valveDisplacement;
            valveBottomPosition_[valveI] += valveDisplacement;
            valvePistonPosition_[valveI] += deltaZ;
        }
    }

    forAll(valves_,valveI)
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

            if (valves_[valveI].curLift() < valveTopTol_)
            {
                topoChanger_[valveLayerID].disable();
            }
            else
            {
                topoChanger_[valveLayerID].enable();
            }
        }

        if(valves_[valveI].isOpen())
        {
            isReallyClosed_[valveI] = false;
        }

    }

    // Changing topology by hand
    autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

    if (pistonLayerID < 0)
    {
        FatalErrorIn("void engineFvMesh::moveAndMorph()")
            << "Piston modifier not found."
            << abort(FatalError);
    }

    // Work array for new points position.
    pointField newPoints = points();

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



// NEW !

    // Reset the position of layered interfaces

/*
    Move the mesh points
    1) Move all the piston points
*/

    bool poppetDeformation = false;

    if(!realDeformation())
    {
        Info << "verticalValves::update()::Layering mode" << endl;

#       include "movePistonPoints.H"
#       include "moveValvePoints.H"
        movePoints(newPoints);
#       include "poppetDeformation.H"
        newPoints = points();
        setVirtualPositions();
    }
    else
    {
        Info << "verticalValves::update()::Deformation mode" << endl;

#       include "moveAllTogether.H"
        movePoints(mSolver.curPoints());
        setVirtualPositions();
        newPoints = points();
#       include "moveTopOfTheValves.H"
        movePoints(newPoints);
#       include "poppetDeformation.H"
        newPoints = points();
        setVirtualPositions();
    }

    forAll(valves(), valveI)
    {
        Info << "Valve Top Position for valve n. " << valveI + 1 << " = " <<
        valveTopPosition_[valveI] << endl;
        Info << "Valve Bottom Position for valve n. " << valveI + 1 << " = " <<
        valveBottomPosition_[valveI] << endl;
        Info << "Valve Piston Position for valve n. " << valveI + 1 << " = " <<
        valvePistonPosition_[valveI] << endl;
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
        pointField oldPointsNew = oldPoints();

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

                mappedOldPointsNew.map
                (
                    oldPointsNew,
                    topoChangeMap3->pointMap()
                );
                pointField newPoints = allPoints();

                movePoints(mappedOldPointsNew);

                resetMotion();
                setV0();
                movePoints(newPoints);
            }
        }
    }

    return true;
}

