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

#include "simpleEngineTopoFvMesh.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "attachDetach.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "tetDecompositionMotionSolver.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleEngineTopoFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        simpleEngineTopoFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::simpleEngineTopoFvMesh::makeLayersLive()
{
    // Enable layering
    forAll (topoChanger_, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanger_[modI]))
        {
            if (debug)
            {
                Info<< "Enabling layer modifier "
                    << topoChanger_[modI].name() << endl;
            }

            topoChanger_[modI].enable();
        }
        else if (isA<slidingInterface>(topoChanger_[modI]))
        {
            if (debug)
            {
                Info<< "Disabling slider modifier "
                    << topoChanger_[modI].name() << endl;
            }

            topoChanger_[modI].disable();
        }
        else if (isA<attachDetach>(topoChanger_[modI]))
        {
            topoChanger_[modI].enable();
        }
        else
        {
            FatalErrorIn("void Foam::simpleEngineTopoFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanger_[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::simpleEngineTopoFvMesh::makeSlidersLive()
{
    // Enable sliding interface
    forAll (topoChanger_, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanger_[modI]))
        {
            if (debug)
            {
                Info<< "Disabling layer modifier "
                    << topoChanger_[modI].name() << endl;
            }

            topoChanger_[modI].disable();
        }
        else if (isA<slidingInterface>(topoChanger_[modI]))
        {
            if (debug)
            {
                Info<< "Enabling slider modifier "
                    << topoChanger_[modI].name() << endl;
            }

            topoChanger_[modI].enable();
        }
        else if (isA<attachDetach>(topoChanger_[modI]))
        {
            topoChanger_[modI].enable();
        }
        else
        {
            FatalErrorIn("void Foam::simpleEngineTopoFvMesh::makeSlidersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanger_[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::simpleEngineTopoFvMesh::prepareValveDetach()
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
                    "void Foam::simpleEngineTopoFvMesh::prepareValveDetach()"
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


bool Foam::simpleEngineTopoFvMesh::attached() const
{
    bool result = false;

    forAll (topoChanger_, modI)
    {
        if (isA<slidingInterface>(topoChanger_[modI]))
        {
            result =
                result
             || refCast<const slidingInterface>(topoChanger_[modI]).attached();
        }
    }

    // Check thal all sliders are in sync (debug only)
    forAll (topoChanger_, modI)
    {
        if (isA<slidingInterface>(topoChanger_[modI]))
        {
            if
            (
                result 
             != refCast<const slidingInterface>(topoChanger_[modI]).attached()
            )
            {
                FatalErrorIn("bool simpleEngineTopoFvMesh::attached() const")
                    << "Slider " << modI << " named "
                    << topoChanger_[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    if (debug)
    {
        if (result)
        {
            Info << "simpleEngineTopoFvMesh is attached" << endl;
        }
        else
        {
            Info << "simpleEngineTopoFvMesh is detached" << endl;
        }
    }

    return result;
}


void Foam::simpleEngineTopoFvMesh::setBoundaryMotion()
{
    // Set the boundary conditions on motion field in order to solve
    // for motion.  Deformation only happens within the cylinder and
    // not in ports - the motion of valve top is set to zero.  Correct
    // using setBoundaryPosition()
    // HJ, 10/Jun/2004  Reconsider
    if (debug)
    {
        Info << "Setting boundary motion" << endl;
    }

    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

    tetPointVectorField& motionU = mSolver.motionU();

    // Set valve velocity
    forAll (valves_, valveI)
    {
        // If valve is present in geometry, set the motion
        if (valves_[valveI].bottomPatchID().active())
        {
            vector valveVel =
                valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

            // Bottom of the valve moves with given velocity
            motionU.boundaryField()[valves_[valveI].bottomPatchID().index()] ==
                valveVel;

            if (debug)
            {
                Info<< "Valve " << valveI << " lift: "
                    << valves_[valveI].curLift()
                    << " velocity: " << valves_[valveI].curVelocity()
                    << endl;
            }
        }

        if (valves_[valveI].poppetPatchID().active())
        {
            // Top of the valve does not move
            motionU.boundaryField()[valves_[valveI].poppetPatchID().index()] ==
                vector::zero;
        }

        if (valves_[valveI].curtainInCylinderPatchID().active())
        {
//             label cicPatchIndex =
//                 valves_[valveI].curtainInCylinderPatchID().index();

//             componentMixedTetPolyPatchVectorField& pf =
//                 refCast<componentMixedTetPolyPatchVectorField>
//                 (
//                     motionU.boundaryField()[cicPatchIndex]
//                 );

//             if (valves_[valveI].isOpen())
//             {
//                 // Get valve coordinate system
//                 const coordinateSystem& vcs = valves_[valveI].cs();
//                 const scalar r = 0.5*valves_[valveI].diameter();

//                 // Get local points in the patch
//                 const pointField& cpGlobal =
//                 motionU.boundaryField()
//                     [cicPatchIndex].patchMesh().localPoints();

//                 pointField cpLocal(vcs.toLocal(cpGlobal));
//                 scalarField mask =
//                     pos
//                     (
//                         cpLocal.component(vector::Z)
//                       + valves_[valveI].curLift()
//                       - valvePosTol_
//                     );

//                 // Calculate motion of valve centre
//                 cpLocal.replace
//                 (
//                     vector::X,
//                     mask*r + (1.0 - mask)*0.0
//                 );

//                 pf.refValue() =
//                 (
//                     mask*(vcs.toGlobal(cpLocal) - cpGlobal)/
//                     engineTime_.deltaT().value()
//                   + (1.0 - mask)*
//                     vector(vcs.axis().x()*valves_[valveI].curVelocity(), 0, 0)
//                 );

//                 pf.valueFraction() =
//                     mask*vector::one + (1.0 - mask)*vector(1, 0, 0);
//             }
//             else
//             {
//                 pf.refValue() = vector::zero;
//                 pf.valueFraction() = vector::one;
//             }
        }
    }

    // Set piston velocity
    if (piston().patchID().active())
    {
        vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

        if (debug)
        {
            Info<< "Piston velocity: " << pistonVel;
        }

        componentMixedTetPolyPatchVectorField& pp =
            refCast<componentMixedTetPolyPatchVectorField>
            (
                motionU.boundaryField()[piston().patchID().index()]
            );

        if (deformation())
        {
            if (debug)
            {
                Info << " deformation"  << endl;
            }

            pp.refValue() = pistonVel;
        }
        else
        {
            if (debug)
            {
                Info << " layering" << endl;
            }

            pp.refValue() = vector::zero;
        }
    }
}


void Foam::simpleEngineTopoFvMesh::setBoundaryPosition()
{
    // Set the boundary position for layer modifiers
    if (debug)
    {
        Info << "Setting boundary position" << endl;
    }

    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

    tetPointVectorField& motionU = mSolver.motionU();

    // Set valve velocity
    forAll (valves_, valveI)
    {
        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

        // If valve is present in geometry, set the motion
        if (valves_[valveI].bottomPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            motionU.boundaryField()[valves_[valveI].bottomPatchID().index()] ==
                valveVel;

            if (debug)
            {
                Info<< "Valve " << valveI << " lift: "
                    << valves_[valveI].curLift()
                    << " velocity: " << valves_[valveI].curVelocity()
                    << endl;
            }

        }

        if (valves_[valveI].poppetPatchID().active())
        {
            // Top of the valve does not move
            motionU.boundaryField()[valves_[valveI].poppetPatchID().index()] ==
                valveVel;
        }
    }

    // Set piston velocity
    if (piston().patchID().active())
    {
        vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

        componentMixedTetPolyPatchVectorField& pp =
            refCast<componentMixedTetPolyPatchVectorField>
            (
                motionU.boundaryField()[piston().patchID().index()]
            );

        if (debug)
        {
            Info<< "Piston velocity: " << pistonVel << endl;
        }

        pp.refValue() = pistonVel;
    }

    motionU.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::simpleEngineTopoFvMesh::simpleEngineTopoFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    engineTime_(refCast<const engineTime>(time())),
    valves_(*this, engineTime_.engineDict().lookup("valves")),
    piston_(*this, engineTime_.engineDict().subDict("piston")),
    msPtr_(motionSolver::New(*this)),
    deformSwitch_(readScalar(engineTime_.engineDict().lookup("deformAngle"))),
    valvePosTol_(readScalar(engineTime_.engineDict().lookup("valvePosTol")))
{
    // Add zones and modifiers if not already there.
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::simpleEngineTopoFvMesh::update()
{
    tetDecompositionMotionSolver& mSolver =
        refCast<tetDecompositionMotionSolver>(msPtr_());

    // Detaching the interface
    if (attached())
    {
        if (debug)
        {
            Info << "Decoupling sliding interfaces" << endl;
        }

        makeSlidersLive();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap1 = topoChanger_.changeMesh();

        if (topoChangeMap1->morphing())
        {
            mSolver.updateMesh(topoChangeMap1());
        }
    }
    else
    {
        if (debug)
        {
            Info << "Sliding interfaces decoupled" << endl;
        }
    }

    // Perform layer action and mesh motion
    makeLayersLive();

    if (debug)
    {
        Info << "Executing layer action" << endl;
    }

    // Find piston mesh modifier
    const label pistonLayerID =
        topoChanger_.findModifierID("pistonLayer");

    if (pistonLayerID < 0)
    {
        FatalErrorIn("void simpleEngineTopoFvMesh::moveAndMorph()")
            << "Piston modifier not found."
            << abort(FatalError);
    }

    if (deformation())
    {
        // Dectivate piston layer
        if (debug)
        {
            Info << "Disabling piston layer (deformation)"<< endl;
        }

        topoChanger_[pistonLayerID].disable();
    }
    else
    {
        // Activate piston layer
        if (debug)
        {
            Info << "Enabling piston layer (deformation)"<< endl;
        }

        topoChanger_[pistonLayerID].enable();
    }

    // Changing topology by hand
    {
        autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

        if (topoChangeMap2->morphing())
        {
            mSolver.updateMesh(topoChangeMap2());

            if (debug)
            {
                Info << "Topology change; executing pre-motion" << endl;
            }

            movePoints(topoChangeMap2->preMotionPoints());
            setV0();
            resetMotion();

        }
    }

    if (deformation())
    {
        if (debug)
        {
            Info << "Mesh deformation mode" << endl;
        }

        setBoundaryMotion();

        // Solve for motion
        mSolver.solve();

        // Dectivate piston layer
        if (debug)
        {
            Info << "Disabling piston layer (topo 2)"<< endl;
        }

        topoChanger_[pistonLayerID].disable();
    }
    else
    {
        if (debug)
        {
            Info << "Piston layering mode" << endl;
        }

        bool tiltedValves = true;

        if (tiltedValves)
        {
            setBoundaryMotion();

            // Solve for motion
            mSolver.solve();
        }

        // Blocking vertical motion
        mSolver.motionU().internalField().replace(vector::Z, 0);

        // Activate piston layer
        if (debug)
        {
            Info << "Enabling piston layer (topo 2)"<< endl;
        }

        topoChanger_[pistonLayerID].enable();
    }

    // Reset the position of layered interfaces
    setBoundaryPosition();

    movePoints(mSolver.curPoints());

    // Attach the interface
    if (debug)
    {
        Info << "Coupling sliding interfaces" << endl;
    }

    makeSlidersLive();
    prepareValveDetach();

    // Changing topology by hand
    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();

        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        if (debug)
        {
            Info << "Moving points post slider attach" << endl;
        }

        if (topoChangeMap3->morphing())
        {
            mSolver.updateMesh(topoChangeMap3());

            if (debug)
            {
                Info << "Moving points post slider attach" << endl;
            }

            pointField newPoints = allPoints();
            pointField mappedOldPointsNew(newPoints.size());

            mappedOldPointsNew.map(oldPointsNew, topoChangeMap3->pointMap());

            // Solve the correct mesh motion to make sure motion fluxes
            // are solved for and not mapped
            movePoints(mappedOldPointsNew);
            resetMotion();
            setV0();
            movePoints(newPoints);
        }
    }

    if (debug)
    {
        Info << "Sliding interfaces coupled: " << attached() << endl;
    }

    return true;
}


// ************************************************************************* //
