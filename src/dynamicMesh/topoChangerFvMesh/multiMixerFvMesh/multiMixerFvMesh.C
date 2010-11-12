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

#include "multiMixerFvMesh.H"
#include "Time.H"
#include "regionSplit.H"
#include "slidingInterface.H"
#include "mapPolyMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMixerFvMesh, 0);
    addToRunTimeSelectionTable(topoChangerFvMesh, multiMixerFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiMixerFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if
    (
        pointZones().size() > 0
     || faceZones().size() > 0
     || cellZones().size() > 0
    )
    {
        Info<< "void multiMixerFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        if (topoChanger_.size() == 0  && useTopoSliding())
        {
            FatalErrorIn
            (
                "void multiMixerFvMesh::addZonesAndModifiers()"
            )   << "Mesh modifiers not read properly"
                << abort(FatalError);
        }

        return;
    }

    Info<< "Time = " << time().timeName() << endl
        << "Adding zones and modifiers to the mesh.  " << rotors_.size()
        << " sliders found" << endl;

    DynamicList<pointZone*> pz(rotors_.size());
    DynamicList<faceZone*> fz(3*rotors_.size());
    DynamicList<cellZone*> cz(rotors_.size());

    // Create region split: mark every cell with its topological region
    regionSplit rs(*this);

    Info << "Adding point, face and cell zones" << endl;
    forAll (rotors_, rotorI)
    {
        rotors_[rotorI].addZones(pz, fz, cz, rs);
    }

    {
        List<pointZone*> pzList;
        pzList.transfer(pz.shrink());

        List<faceZone*> fzList;
        fzList.transfer(fz.shrink());

        List<cellZone*> czList;
        czList.transfer(cz.shrink());

        addZones(pzList, fzList, czList);
    }

    if (useTopoSliding())
    {
        topoChanger_.setSize(rotors_.size());
        label nextI = 0;

        forAll (rotors_, rotorI)
        {
            rotors_[rotorI].addModifiers(topoChanger_, nextI);
        }

        Info<< "Adding topology modifiers.  nModifiers = " << nextI << endl;

        // Resize and set topo changer
        topoChanger_.setSize(nextI);
        topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
        topoChanger_.write();
    }

    // Write mesh with new zones
    write();
}


bool Foam::multiMixerFvMesh::useTopoSliding() const
{
    bool result = false;

    forAll (rotors_, rotorI)
    {
        result = (result || rotors_[rotorI].useTopoSliding());
    }

    return result;
}


void Foam::multiMixerFvMesh::checkRotors() const
{
    // Check if all sliding interfaces are in the same state
    bool rotorState = attached();

    forAll (topoChanger_, topoI)
    {
        if (isA<slidingInterface>(topoChanger_[topoI]))
        {
            const slidingInterface& slider =
                refCast<const slidingInterface>(topoChanger_[topoI]);

            if (rotorState != slider.attached())
            {
                FatalErrorIn
                (
                    "void multiMixerFvMesh::checkRotors() const"
                )   << "Some sliders are attached and some are not.  "
                    << "Out of sync sliders cannot be handled."
                    << abort(FatalError);
            }
        }
    }
}


bool Foam::multiMixerFvMesh::attached() const
{
    // Check if interfaces are attached AND if they are in sync
    bool attached = false;

    forAll (topoChanger_, topoI)
    {
        if (isA<slidingInterface>(topoChanger_[topoI]))
        {
            attached = attached
                || refCast<const slidingInterface>(topoChanger_[0]).attached();
        }
    }

    return attached;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::multiMixerFvMesh::multiMixerFvMesh
(
    const IOobject& io
)
:
    topoChangerFvMesh(io),
    dict_
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
    ),
    rotors_(),
    attachDetach_(dict_.lookupOrDefault<bool>("attachDetach", true))
{
    // Read rotors from the dictionary
    PtrList<entry> rotorEntries(dict_.lookup("rotors"));
    rotors_.setSize(rotorEntries.size());

    forAll (rotorEntries, rotorI)
    {
        rotors_.set
        (
            rotorI,
            new mixerRotor
            (
                rotorEntries[rotorI].keyword(),
                *this,
                rotorEntries[rotorI].dict()
            )
        );
    }

    addZonesAndModifiers();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiMixerFvMesh::~multiMixerFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::multiMixerFvMesh::update()
{
    if (useTopoSliding())
    {
        // Detaching the interface
        if (attached() && attachDetach_)
        {
            Info << "Detaching rotors" << endl;
            autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

            forAll (rotors_, rotorI)
            {
                rotors_[rotorI].updateTopology();
            }
        }
    }

    // Accumulate point motion
    vectorField pointMotion(allPoints().size(), vector::zero);

    forAll (rotors_, rotorI)
    {
        pointMotion += rotors_[rotorI].pointMotion();
    }

    // Save old points
    pointField oldPointsNew = allPoints();

    // Move points
    movePoints(allPoints() + pointMotion);

    if (useTopoSliding())
    {
        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();
        bool morphing = topoChangeMap.valid();

        if (morphing)
        {
            Info << "Attaching rotors" << endl;

            forAll (rotors_, rotorI)
            {
                rotors_[rotorI].updateTopology();
            }

            // Move the sliding interface points to correct position
            pointField mappedOldPointsNew(allPoints().size());
            mappedOldPointsNew.map(oldPointsNew, topoChangeMap->pointMap());

            movePoints(mappedOldPointsNew);
            resetMotion();
            setV0();

            // Move the sliding interface points to correct position
            movePoints(topoChangeMap->preMotionPoints());
        }
    }

    return true;
}


// ************************************************************************* //
