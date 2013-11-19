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

#include "linearValveLayersFvMesh.H"
#include "Time.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "pointField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearValveLayersFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        linearValveLayersFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::linearValveLayersFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "void linearValveLayersFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0)
        {
            FatalErrorIn
            (
                "void linearValveLayersFvMesh::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    List<pointZone*> pz(1);
    List<faceZone*> fz(4);
    label nFz = 0;
    List<cellZone*> cz(0);

    // Add a topology modifier
    topoChanger_.setSize(2);
    label nTc = 0;

    // Add an empty zone for cut points

    pz[0] = new pointZone
    (
        "cutPointZone",
        labelList(0),
        0,
        pointZones()
    );


    // Do face zones for slider

    // Inner slider
    const word innerSliderName(motionDict_.subDict("slider").lookup("inside"));
    const polyPatch& innerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(innerSliderName)];

    // Outer slider
    const word outerSliderName
    (
        motionDict_.subDict("slider").lookup("outside")
    );

    const polyPatch& outerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(outerSliderName)];

    if (!innerSlider.empty() && !outerSlider.empty())
    {
        Pout<< "Adding sliding interface between patches "
            << innerSliderName << " and " << outerSliderName << endl;
        labelList isf(innerSlider.size());

        forAll (isf, i)
        {
            isf[i] = innerSlider.start() + i;
        }

        fz[0] = new faceZone
        (
            innerSliderName + "Zone",
            isf,
            boolList(innerSlider.size(), false),
            0,
            faceZones()
        );

        labelList osf(outerSlider.size());

        forAll (osf, i)
        {
            osf[i] = outerSlider.start() + i;
        }

        fz[1] = new faceZone
        (
            outerSliderName + "Zone",
            osf,
            boolList(outerSlider.size(), false),
            1,
            faceZones()
        );

        // Add empty zone for cut faces
        fz[2] = new faceZone
        (
            "cutFaceZone",
            labelList(0),
            boolList(0, false),
            2,
            faceZones()
        );

        // Set the number of face zones already used
        nFz = 3;
    }

    // Add face zone for layer addition
    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("patch")
    );

    const polyPatch& layerPatch =
        boundaryMesh()[boundaryMesh().findPatchID(layerPatchName)];

    labelList lpf(layerPatch.size());

    forAll (lpf, i)
    {
        lpf[i] = layerPatch.start() + i;
    }

    fz[nFz] = new faceZone
    (
        "valveLayerZone",
        lpf,
        boolList(layerPatch.size(), true),
        0,
        faceZones()
    );
    nFz++;
    Pout<< "nFz = " << nFz << endl;
    // Resize the number of live face zones
    fz.setSize(nFz);


    Info << "Adding point and face zones" << endl;
    addZones(pz, fz, cz);

    if (!innerSlider.empty() && !outerSlider.empty())
    {
        // Set the topo changer for sliding interface
        topoChanger_.set
        (
            0,
            new slidingInterface
            (
                "valveSlider",
                0,
                topoChanger_,
                outerSliderName + "Zone",
                innerSliderName + "Zone",
                "cutPointZone",
                "cutFaceZone",
                outerSliderName,
                innerSliderName,
                slidingInterface::INTEGRAL,   // Edge matching algorithm
                true,                         // Attach-detach action
                intersection::VISIBLE         // Projection algorithm
            )
        );

        // Record one added topo modifyer
        nTc = 1;
    }

    // Add the layer addition-removal topo modifyer
    topoChanger_.set
    (
        nTc,
        new layerAdditionRemoval
        (
            "valveLayer",
            nTc,
            topoChanger_,
            "valveLayerZone",
            readScalar
            (
                motionDict_.subDict("layer").lookup("minThickness")
            ),
            readScalar
            (
                motionDict_.subDict("layer").lookup("maxThickness")
            )
        )
    );
    nTc++;

    // Reset the size of mesh modifiers
    topoChanger_.setSize(nTc);
    Pout<< nTc << " topology modifiers on processor" << endl;
    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}


void Foam::linearValveLayersFvMesh::makeLayersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable layering
    forAll (topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else
        {
            FatalErrorIn("void linearValveLayersFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


void Foam::linearValveLayersFvMesh::makeSlidersLive()
{
    const polyTopoChanger& topoChanges = topoChanger_;

    // Enable sliding interface
    forAll (topoChanges, modI)
    {
        if (isA<layerAdditionRemoval>(topoChanges[modI]))
        {
            topoChanges[modI].disable();
        }
        else if (isA<slidingInterface>(topoChanges[modI]))
        {
            topoChanges[modI].enable();
        }
        else
        {
            FatalErrorIn("void linearValveLayersFvMesh::makeLayersLive()")
                << "Don't know what to do with mesh modifier "
                << modI << " of type " << topoChanges[modI].type()
                << abort(FatalError);
        }
    }
}


bool Foam::linearValveLayersFvMesh::attached() const
{
    const polyTopoChanger& topoChanges = topoChanger_;

    bool result = false;

    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            result =
                result
             || refCast<const slidingInterface>(topoChanges[modI]).attached();
        }
    }

    // Check that all sliders are in sync (debug only)
    forAll (topoChanges, modI)
    {
        if (isA<slidingInterface>(topoChanges[modI]))
        {
            if
            (
                result
             != refCast<const slidingInterface>(topoChanges[modI]).attached()
            )
            {
                FatalErrorIn("bool linearValveLayersFvMesh::attached() const")
                    << "Slider " << modI << " named "
                    << topoChanges[modI].name()
                    << " out of sync: Should be" << result
                    << abort(FatalError);
            }
        }
    }

    // Sync across processors
    reduce(result, orOp<bool>());

    return result;
}


Foam::tmp<Foam::pointField> Foam::linearValveLayersFvMesh::newPoints() const
{
    tmp<pointField> tnewPoints
    (
        new pointField(points())
    );

    pointField& np = tnewPoints();

    const word layerPatchName
    (
        motionDict_.subDict("layer").lookup("patch")
    );

    const polyPatch& layerPatch =
        boundaryMesh()[boundaryMesh().findPatchID(layerPatchName)];

    const labelList& patchPoints = layerPatch.meshPoints();

    const vector vel
    (
        motionDict_.lookup("pistonVelocity")
    );

    forAll (patchPoints, ppI)
    {
        np[patchPoints[ppI]] += vel*time().deltaT().value();
    }

    return tnewPoints;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::linearValveLayersFvMesh::linearValveLayersFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    motionDict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    )
{
    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearValveLayersFvMesh::~linearValveLayersFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::linearValveLayersFvMesh::update()
{
    // Detaching the interface
    if (attached())
    {
        Info<< "Decoupling sliding interfaces" << endl;
        makeSlidersLive();

        // Changing topology by hand
        autoPtr<mapPolyMesh> topoChangeMap1 = topoChanger_.changeMesh();

        bool localMorphing1 = topoChangeMap1->morphing();

        if (localMorphing1)
        {
            Info << "Topology change; executing pre-motion after "
                << "sliding detach" << endl;
            movePoints(topoChangeMap1->preMotionPoints());
        }
        else
        {
            pointField newPoints = allPoints();

            // Dummy motion
            movePoints(newPoints);
        }

        Info<< "sliding interfaces successfully decoupled!!!" << endl;
    }
    else
    {
        Info<< "Sliding interfaces decoupled" << endl;
    }

    // Perform layer action and mesh motion
    makeLayersLive();

    // Changing topology by hand
    {
        autoPtr<mapPolyMesh> topoChangeMap2 = topoChanger_.changeMesh();

        bool localMorphing2 = topoChangeMap2->morphing();
        bool globalMorphing2 = localMorphing2;

        reduce(globalMorphing2, orOp<bool>());

        // Work array for new points position.
        pointField newPoints = allPoints();

        if (globalMorphing2)
        {
            Info<< "Topology change; executing pre-motion after "
                << "dynamic layering" << endl;

            if (localMorphing2)
            {
                Info << "Topology change; executing pre-motion" << endl;
                movePoints(topoChangeMap2->preMotionPoints());
            }
            else
            {
                movePoints(newPoints);
            }
        }
    }

    // Move points to change layer thickness
    movePoints(newPoints());


    // Attach the interface
    Info << "Coupling sliding interfaces" << endl;
    makeSlidersLive();

    // Changing topology by hand
    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();

        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        bool localMorphing3 = topoChangeMap3->morphing();
        bool globalMorphing3 = localMorphing3;

        reduce(globalMorphing3, orOp<bool>());

        pointField newPoints = allPoints();

        if (globalMorphing3)
        {
            if (localMorphing3)
            {
                pointField mappedOldPointsNew(newPoints.size());

                mappedOldPointsNew.map
                (
                    oldPointsNew,
                    topoChangeMap3->pointMap()
                );

                // Solve the correct mesh motion to make sure motion fluxes
                // are solved for and not mapped
                movePoints(mappedOldPointsNew);
            }
            else
            {
                movePoints(newPoints);
            }
        }

        // Reset motion
        resetMotion();
        setV0();
        movePoints(newPoints);
    }

    return true;
}


// ************************************************************************* //
