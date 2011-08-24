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

#include "movingBodyTopoFvMesh.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "layerAdditionRemoval.H"
#include "volMesh.H"
#include "transformField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(movingBodyTopoFvMesh, 0);

    addToRunTimeSelectionTable
    (
        topoChangerFvMesh,
        movingBodyTopoFvMesh,
        IOobject
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::movingBodyTopoFvMesh::calcMotionMask() const
{
    Info<< "Updating vertex markup" << endl;

    tmp<scalarField> tvertexMarkup(new scalarField(allPoints().size(), 0));
    scalarField& vertexMarkup = tvertexMarkup();

    cellZoneID movingCellsID(movingCellsName_, cellZones());

    // In order to do a correct update on a mask on processor boundaries,
    // Detection of moving cells should use patchNeighbourField for
    // processor (not coupled!) boundaries.  This is done by expanding
    // a moving cell set into a field and making sure that processor patch
    // points move in sync.  Not done at the moment, probably best to do
    // using parallel update of pointFields.  HJ, 19/Feb/2011

    // If moving cells are found, perform mark-up
    if (movingCellsID.active())
    {
        // Get cell-point addressing
        const labelListList& cp = cellPoints();

        // Get labels of all moving cells
        const labelList& movingCells = cellZones()[movingCellsID.index()];

        forAll (movingCells, cellI)
        {
            const labelList& curCp = cp[movingCells[cellI]];

            forAll (curCp, pointI)
            {
                vertexMarkup[curCp[pointI]] = 1;
            }
        }
    }

    faceZoneID frontFacesID(frontFacesName_, faceZones());

    if (frontFacesID.active())
    {
        const faceZone& frontFaces = faceZones()[frontFacesID.index()];

        const labelList& mp = frontFaces().meshPoints();

        forAll (mp, mpI)
        {
            vertexMarkup[mp[mpI]] = 1;
        }
    }

    faceZoneID backFacesID(backFacesName_, faceZones());

    if (backFacesID.active())
    {
        const faceZone& backFaces = faceZones()[backFacesID.index()];

        const labelList& mp = backFaces().meshPoints();

        forAll (mp, mpI)
        {
            vertexMarkup[mp[mpI]] = 1;
        }
    }

    return tvertexMarkup;
}


void Foam::movingBodyTopoFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action

    if (topoChanger_.size() > 0)
    {
        Info<< "void movingBodyTopoFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }

    // Add layer addition/removal interfaces
    topoChanger_.setSize(2);
    label nMods = 0;


    faceZoneID frontFacesID(frontFacesName_, faceZones());
    faceZoneID backFacesID(backFacesName_, faceZones());

    if (frontFacesID.active())
    {
        const faceZone& frontFaces = faceZones()[frontFacesID.index()];

        if (!frontFaces.empty())
        {
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    frontFacesName_ + "Layer",
                    nMods,
                    topoChanger_,
                    frontFacesName_,
                    readScalar
                    (
                        dict_.subDict("front").lookup("minThickness")
                    ),
                    readScalar
                    (
                        dict_.subDict("front").lookup("maxThickness")
                    )
                )
            );

            nMods++;
        }
    }

    if (backFacesID.active())
    {
        const faceZone& backFaces = faceZones()[backFacesID.index()];

        if (!backFaces.empty())
        {
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    backFacesName_ + "Layer",
                    nMods,
                    topoChanger_,
                    backFacesName_,
                    readScalar
                    (
                        dict_.subDict("back").lookup("minThickness")
                    ),
                    readScalar
                    (
                        dict_.subDict("back").lookup("maxThickness")
                    )
                )
            );

            nMods++;
        }
    }

    topoChanger_.setSize(nMods);

    reduce(nMods, sumOp<label>());

    Info << "Adding " << nMods << " mesh modifiers" << endl;

    // Write mesh and modifiers
    topoChanger_.write();

    // No need to write the mesh - only modifiers are added.
    // HJ, 18/Feb/2011
//     write();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::movingBodyTopoFvMesh::movingBodyTopoFvMesh(const IOobject& io)
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
    movingCellsName_(dict_.lookup("movingCells")),
    frontFacesName_(dict_.lookup("frontFaces")),
    backFacesName_(dict_.lookup("backFaces")),
    SBMFPtr_(solidBodyMotionFunction::New(dict_, time())),
    motionMask_()
{
    addZonesAndModifiers();
    motionMask_ = calcMotionMask();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::movingBodyTopoFvMesh::~movingBodyTopoFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::movingBodyTopoFvMesh::update()
{
    // Store points to recreate mesh motion
    pointField oldPointsNew = allPoints();
    pointField newPoints = allPoints();

    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    bool localMeshChanged = topoChangeMap->morphing();
    bool globalMeshChanged = localMeshChanged;
    reduce(globalMeshChanged, orOp<bool>());

    if (globalMeshChanged)
    {
        Pout<< "Topology change. Calculating motion point mask" << endl;
        motionMask_ = calcMotionMask();
    }

    if (localMeshChanged)
    {
//         // Map old points onto the new mesh
//         pointField mappedOldPointsNew(allPoints().size());
//         mappedOldPointsNew.map(oldPointsNew, topoChangeMap->pointMap());

//         movePoints(mappedOldPointsNew);
//         resetMotion();
//         setV0();

        // Get new points from preMotion
        newPoints = topoChangeMap().preMotionPoints();
    }
//     else
//     {
//         // No change, use old points
//         movePoints(oldPointsNew);
//         resetMotion();
//         setV0();
//     }

    // Calculate new points using a velocity transformation
    newPoints += motionMask_*
        transform(SBMFPtr_().velocity(), newPoints)*time().deltaT().value();

    Info << "Executing mesh motion" << endl;
    movePoints(newPoints);

    return localMeshChanged;
}


// ************************************************************************* //
