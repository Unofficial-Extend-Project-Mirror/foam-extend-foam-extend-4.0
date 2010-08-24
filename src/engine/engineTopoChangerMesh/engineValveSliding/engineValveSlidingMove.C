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

#include "engineValveSliding.H"
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

void Foam::engineValveSliding::makeLayersLive()
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

void Foam::engineValveSliding::makeSlidersLive()
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

bool Foam::engineValveSliding::attached() const
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


void Foam::engineValveSliding::valveDetach()
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

void Foam::engineValveSliding::valveAttach()
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


void Foam::engineValveSliding::prepareValveDetach()
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


bool Foam::engineValveSliding::update()
{

    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

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

    scalar deltaZ = engTime().pistonDisplacement().value();

    // Piston Layering
    makeLayersLive();

    {
        forAll(valves_, valveI)
        {
            scalar valveDisplacement = valves_[valveI].curVelocity()*
                valves_[valveI].cs().axis().z()*engTime().deltaT().value();

            Info << "valveDisplacement = " << valveDisplacement << nl
                << "valve velocity =" << valves_[valveI].curVelocity() << nl
                << "valve lift =" << valves_[valveI].curLift() << nl
                << "valve deformation lift ="
                << valves_[valveI].deformationLift() << endl;
        }

    }

    forAll(valves_,valveI)
    {
        if(valves_[valveI].curLift() >= valves_[valveI].deformationLift())
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


#        include "setMotionBoundaryConditionEngineValveSliding.H"

        DynamicList<label> constrainedPoints(mSolver.curPoints()().size()/100);
        DynamicList<vector> constrainedVelocity
        (
            mSolver.curPoints()().size()/100
        );

#       include "setEngineValveSlidingConstraint.H"

        forAll (constrainedPoints, i)
        {
            mSolver.setConstraint
            (
                constrainedPoints[i],
                constrainedVelocity[i]
            );
        }

//        Info << mSolver.motionU() << endl;

        Info << "pre solve" << endl;

        mSolver.solve();

        newPoints = mSolver.curPoints();
        movePoints(newPoints);

        Info << "max mesh phi, 1 = " << max(phi()) << endl;
        Info << "min mesh phi, 1 = " << min(phi()) << endl;

        setVirtualPositions();
        mSolver.clearConstraints();



//      layering

#       include "moveValvePointsEngineValveSliding.H"
        movePoints(newPoints);
//        newPoints = points();
        newPoints = allPoints();
        setVirtualPositions();

        Info << "max mesh phi, 2 = " << max(phi()) << endl;
        Info << "min mesh phi, 2 = " << min(phi()) << endl;

    }

    Info << ", deckHeight = " << deckHeight()
        << ", pistonPosition = " << pistonPosition()    << endl;

    pistonPosition() += deltaZ;

    Info << "max V()-V0() before = " << max(mag(V() - V0())) << endl;
    Info << "max V() before = " << max(mag(V())) << endl;

    Info<< "BEFORE mag(points - oldPoints)/deltaT = "
        << max(mag(points() - oldPoints()))/engTime().deltaT().value() << endl;

/*
    forAll(points(), pp)
    {
        Info<< "OLD " <<  "point " << pp << " " << points()[pp]
            << " old point = " << oldPoints()[pp] << mag(points()[pp] -
        oldPoints()[pp])/engTime().deltaT().value() << endl;
    }

*/

//*/ //Tommaso, 23/5/2008

    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();
//        pointField oldPointsNew = oldPoints();

        // Attach the interface
        Info << "Coupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap4 = topoChanger_.changeMesh();

        // Changing topology by hand
        if(topoChangeMap4->morphing())
        {
            mSolver.updateMesh(topoChangeMap4());

            if (topoChangeMap4->hasMotionPoints())
            {
//                 Info<< "Topology change; executing pre-motion "
//                     << "after sliding attach" << endl;

//                 Info<< "topoChangeMap4->preMotionPoints().size() = "
//                     << topoChangeMap4->preMotionPoints().size() << endl;
//                 Info << "allPoints.size() = " << allPoints().size() << endl;
//                 Info << "points.size() = " << points().size() << endl;

                movePoints(topoChangeMap4->preMotionPoints());
//            newPoints = topoChangeMap4->preMotionPoints();
//            newPoints = allPoints();
                newPoints = allPoints();
            }
        }

        surfaceScalarField correctedMeshPhi = phi();


        Info << "max V()-V0() after = " << max(mag(V() - V0())) << endl;
        Info << "max V() before = " << max(mag(V())) << endl;

        {
            Info << "AFTER mag(points - oldPoints)/deltaT = "
                << max(mag(points() - oldPoints()))/engTime().deltaT().value()
                << endl;
            Info << "max V()-V0() after movePoints?!? = "
                << max(mag(V() - V0())) << endl;

            pointField mappedOldPointsNew(allPoints().size());

            mappedOldPointsNew.map(oldPointsNew, topoChangeMap4->pointMap());

            forAll(valves(), valveI)
            {
                const labelList& cutPointsAddressing =
                    pointZones()
                    [
                        pointZones().findZoneID
                        (
                            "cutPointsV"
                            + Foam::name(valveI + 1)
                        )
                    ];

                forAll(cutPointsAddressing, i)
                {
                    mappedOldPointsNew[cutPointsAddressing[i]] =
                        newPoints[cutPointsAddressing[i]];
                }
            }

            pointField newPoints = allPoints();

            movePoints(mappedOldPointsNew);

            resetMotion();
            setV0();
            movePoints(newPoints);

            Info<< "max mesh phi, sliding mesh attached = "
                << max(mag(phi().internalField())) << endl;
            Info<< "max mesh phi, sliding mesh attached = "
                << min(mag(phi().internalField())) << endl;

            Info<< "AFTER 2 mag(points - oldPoints)/deltaT = "
                << max(mag(points() - oldPoints()))/engTime().deltaT().value()
                << endl;
            Info << "AFTER 2 max V()-V0() after = "
                << max(mag(V() - V0())) << endl;

        }
    }

    return true;
}


