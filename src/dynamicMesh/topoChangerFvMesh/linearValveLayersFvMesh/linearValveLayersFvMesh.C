/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "linearValveLayersFvMesh.H"
#include "foamTime.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "pointField.H"
#include "mapPolyMesh.H"
#include "polyTopoChange.H"
#include "pointZone.H"
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
    // Inner slider
    const word innerSliderName(motionDict_.subDict("slider").lookup("inside"));

    // Outer slider
    const word outerSliderName
    (
        motionDict_.subDict("slider").lookup("outside")
    );

    bool initialised = false;

    // Check if zones and modifiers for motion action are present
    label insideZoneID = faceZones().findZoneID(innerSliderName + "Zone");
    label outsideZoneID = faceZones().findZoneID(outerSliderName + "Zone");

    if
    (
        insideZoneID > -1
     || outsideZoneID > -1
    )
    {
        // Zones found.  Check topo changer

        if (topoChanger_.empty())
        {
            FatalErrorIn
            (
                "void linearValveLayersFvMesh::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        initialised = true;
    }

    // Check if slider has been initialised on any of the processors
    reduce(initialised, orOp<bool>());

    if (initialised)
    {
        InfoIn("void linearValveLayersFvMesh::addZonesAndModifiers()")
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    // Add zones and modifiers for motion action
    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh" << endl;

    // Add zones
    label nPz = 0;
    label nFz = 0;
    label nCz = 0;
    List<pointZone*> pz(pointZones().size() + 1);
    List<faceZone*> fz(faceZones().size() + 4);
    List<cellZone*> cz(cellZones().size());

    // Add a topology modifier
    topoChanger_.setSize(2);
    label nTc = 0;

    // Copy existing point zones
    forAll (pointZones(), zoneI)
    {
        pz[nPz] = pointZones()[zoneI].clone(pointZones()).ptr();
        nPz++;
    }

    // Copy existing face zones
    forAll (faceZones(), zoneI)
    {
        fz[nFz] = faceZones()[zoneI].clone(faceZones()).ptr();
        nFz++;
    }

    // Copy existing cell zones
    forAll (cellZones(), zoneI)
    {
        cz[nCz] = cellZones()[zoneI].clone(cellZones()).ptr();
        nCz++;
    }



    // Do face zones for slider

    // Inner slider
    const polyPatch& innerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(innerSliderName)];

    // Outer slider
    const polyPatch& outerSlider =
        boundaryMesh()[boundaryMesh().findPatchID(outerSliderName)];

    if (!innerSlider.empty() && !outerSlider.empty())
    {
        Pout<< "Adding sliding interface between patches "
            << innerSliderName << " and " << outerSliderName << endl;


        // Add an empty zone for cut points
        pz[nPz] = new pointZone
        (
            "cutPointZone",
            labelList(0),
            nPz,
            pointZones()
        );
        nPz++;

        labelList isf(innerSlider.size());

        forAll (isf, i)
        {
            isf[i] = innerSlider.start() + i;
        }

        fz[nFz] = new faceZone
        (
            innerSliderName + "Zone",
            isf,
            boolList(innerSlider.size(), false),
            nFz,
            faceZones()
        );
        nFz++;

        labelList osf(outerSlider.size());

        forAll (osf, i)
        {
            osf[i] = outerSlider.start() + i;
        }

        fz[nFz] = new faceZone
        (
            outerSliderName + "Zone",
            osf,
            boolList(outerSlider.size(), false),
            nFz,
            faceZones()
        );
        nFz++;

        // Add empty zone for cut faces
        fz[nFz] = new faceZone
        (
            "cutFaceZone",
            labelList(0),
            boolList(0, false),
            nFz,
            faceZones()
        );
        nFz++;
    }

    // Add face zone for layer addition.  This is present on all processors
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
        nFz,
        faceZones()
    );
    nFz++;

    // Resize the number of live zones
    pz.setSize(nPz);
    fz.setSize(nFz);
    // Cell zones remain unchanged

    Info << "Adding point and face zones" << endl;
    removeZones();
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

    // Update the mesh for changes in zones.  This needs to
    // happen after write, because mesh instance will be changed
    // HJ and OP, 20/Nov/2013
    syncUpdateMesh();
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


Foam::tmp<Foam::pointField>
Foam::linearValveLayersFvMesh::newLayerPoints() const
{
    tmp<pointField> tnewLayerPoints
    (
        new pointField(allPoints())
    );

    pointField& np = tnewLayerPoints();

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

    return tnewLayerPoints;
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

        // Note: Since we are detaching, global morphing is always true
        // HJ, 7/Mar/2011

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
            // Note: using setOldPoints instead of movePoints.
            // HJ, 23/Aug/2015
            setOldPoints(topoChangeMap2->preMotionPoints());
            newPoints = topoChangeMap2->preMotionPoints();
        }
        else
        {
            // Note: using setOldPoints instead of movePoints.
            // HJ, 23/Aug/2015
            setOldPoints(newPoints);
        }

        setV0();
        resetMotion();
    }

    // Move points to change layer thickness
    movePoints(newLayerPoints());

    // Changing topology by hand
    {
        // Grab old points to correct the motion
        pointField oldPointsNew = oldAllPoints();

        // Attach the interface
        Info << "Coupling sliding interfaces" << endl;
        makeSlidersLive();

        autoPtr<mapPolyMesh> topoChangeMap3 = topoChanger_.changeMesh();

        bool localMorphing3 = topoChangeMap3->morphing();
        bool globalMorphing3 = localMorphing3;

        reduce(globalMorphing3, orOp<bool>());

        if (globalMorphing3)
        {
            Info<< "Topology change; executing pre-motion after "
                << "sliding attach" << endl;

            // Grab points
            newPoints = allPoints();

            if (localMorphing3)
            {
                // If there is layering, pick up correct points
                if (topoChangeMap3->hasMotionPoints())
                {
                    newPoints = topoChangeMap3->preMotionPoints();
                }

                pointField mappedOldPointsNew(newPoints.size());

                mappedOldPointsNew.map
                (
                    oldPointsNew,
                    topoChangeMap3->pointMap()
                );

                // Solve the correct mesh motion to make sure motion fluxes
                // are solved for and not mapped
                // Note: using setOldPoints instead of movePoints.
                // HJ, 23/Aug/2015
                setOldPoints(mappedOldPointsNew);

                resetMotion();
                setV0();

                // Set new point motion
                movePoints(newPoints);
            }
            else
            {
                // No local topological change.  Execute double motion for
                // sync with topological changes
                // Note: using setOldPoints instead of movePoints.
                // HJ, 23/Aug/2015
                setOldPoints(oldPointsNew);

                resetMotion();
                setV0();

                // Set new point motion
                movePoints(newPoints);
            }
        }
    }

    return true;
}


// ************************************************************************* //
