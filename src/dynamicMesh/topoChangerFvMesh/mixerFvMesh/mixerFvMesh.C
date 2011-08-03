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

#include "mixerFvMesh.H"
#include "Time.H"
#include "regionSplit.H"
#include "slidingInterface.H"
#include "mapPolyMesh.H"
#include "volMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixerFvMesh, 0);
    addToRunTimeSelectionTable(topoChangerFvMesh, mixerFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mixerFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action
    if (cellZones().size() > 0)
    {
        Info<< "void mixerFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        // Check definition of the modifier
        if
        (
            pointZones().size() > 0
         || faceZones().size() > 0
        )
        {
            if (topoChanger_.size() == 0)
            {
                FatalErrorIn
                (
                    "void mixerFvMesh::addZonesAndModifiers()"
                )   << "Mesh modifiers not read properly.  "
                    << "pointZones = " <<  pointZones().size()
                    << " faceZones = " << faceZones().size()
                    << abort(FatalError);
            }
        }

        return;
    }

    // Find patches and check sizes

    // Moving slider
    const word movingSliderName(dict_.subDict("slider").lookup("moving"));
    label movingSliderIndex = boundaryMesh().findPatchID(movingSliderName);

    if (movingSliderIndex < 0)
    {
        FatalErrorIn("void mixerFvMesh::addZonesAndModifiers() const")
            << "Moving slider patch not found in boundary"
            << abort(FatalError);
    }

    const word staticSliderName(dict_.subDict("slider").lookup("static"));
    label staticSliderIndex = boundaryMesh().findPatchID(staticSliderName);

    if (staticSliderIndex < 0)
    {
        FatalErrorIn("void mixerFvMesh::addZonesAndModifiers() const")
            << "Static slider patch not found in boundary"
            << abort(FatalError);

    }

    const polyPatch& movingSlider = boundaryMesh()[movingSliderIndex];
    const polyPatch& staticSlider = boundaryMesh()[staticSliderIndex];

    List<pointZone*> pz;
    List<faceZone*> fz;

    bool addSlider = false;

    if (movingSlider.size() > 0 && staticSlider.size() > 0)
    {
        addSlider = true;

        pz.setSize(1);
        fz.setSize(3);

        Info<< "Time = " << time().timeName() << endl
            << "Adding zones and modifiers to the mesh" << endl;

        // Add zones

        // Add an empty zone for cut points

        pz[0] = new pointZone
        (
            "cutPointZone",
            labelList(0),
            0,
            pointZones()
        );


        // Do face zones for slider

        // Moving slider
        labelList isf(movingSlider.size());

        forAll (isf, i)
        {
            isf[i] = movingSlider.start() + i;
        }

        fz[0] = new faceZone
        (
            movingSliderName + "Zone",
            isf,
            boolList(movingSlider.size(), false),
            0,
            faceZones()
        );

        // Static slider
        labelList osf(staticSlider.size());

        forAll (osf, i)
        {
            osf[i] = staticSlider.start() + i;
        }

        fz[1] = new faceZone
        (
            staticSliderName + "Zone",
            osf,
            boolList(staticSlider.size(), false),
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
    }

    // Add cell zone even if slider does not exist

    List<cellZone*> cz(1);

    // Mark every cell with its topological region
    regionSplit rs(*this);

    // Get the region of the cell containing the origin.
    label originRegion = rs[findNearestCell(rotatingRegionMarker_)];

    labelList movingCells(nCells());
    label nMovingCells = 0;

    forAll(rs, cellI)
    {
        if (rs[cellI] == originRegion)
        {
            movingCells[nMovingCells] = cellI;
            nMovingCells++;
        }
    }

    movingCells.setSize(nMovingCells);
    Info << "Number of cells in the moving region: " << nMovingCells << endl;

    cz[0] = new cellZone
    (
        "movingCells",
        movingCells,
        0,
        cellZones()
    );

    Info << "Adding point, face and cell zones" << endl;
    addZones(pz, fz, cz);

    if (addSlider)
    {
        // Add a topology modifier
        Info << "Adding topology modifiers" << endl;
        topoChanger_.setSize(1);

        topoChanger_.set
        (
            0,
            new slidingInterface
            (
                "mixerSlider",
                0,
                topoChanger_,
                staticSliderName + "Zone",
                movingSliderName + "Zone",
                "cutPointZone",
                "cutFaceZone",
                staticSliderName,
                movingSliderName,
                slidingInterface::INTEGRAL,   // Edge matching algorithm
                attachDetach_,                // Attach-detach action
                intersection::VISIBLE         // Projection algorithm
            )
        );
    }

    // Write mesh and modifiers
    topoChanger_.writeOpt() = IOobject::AUTO_WRITE;
    topoChanger_.write();
    write();
}


bool Foam::mixerFvMesh::attached() const
{
    bool result = false;

    if (topoChanger_.size() > 0)
    {
        result = refCast<const slidingInterface>(topoChanger_[0]).attached();
    }

    reduce(result, orOp<bool>());

    return result;
}


void Foam::mixerFvMesh::calcMovingMask() const
{
    if (debug)
    {
        Info<< "void mixerFvMesh::calcMovingMask() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsMaskPtr_)
    {
        FatalErrorIn("void mixerFvMesh::calcMovingMask() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = new scalarField(allPoints().size(), 0);
    scalarField& movingPointsMask = *movingPointsMaskPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    const labelList& cellAddr =
        cellZones()[cellZones().findZoneID("movingCells")];

    forAll (cellAddr, cellI)
    {
        const cell& curCell = c[cellAddr[cellI]];

        forAll (curCell, faceI)
        {
            // Mark all the points as moving
            const face& curFace = f[curCell[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }
    }

    if (topoChanger_.size() > 0)
    {
        // Topo changer present.  Use zones for motion
        const word movingSliderZoneName
        (
            word(dict_.subDict("slider").lookup("moving")) + "Zone"
        );

        const labelList& movingSliderAddr =
            faceZones()[faceZones().findZoneID(movingSliderZoneName)];

        forAll (movingSliderAddr, faceI)
        {
            const face& curFace = f[movingSliderAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }

        const word staticSliderZoneName
        (
            word(dict_.subDict("slider").lookup("static")) + "Zone"
        );

        const labelList& staticSliderAddr =
            faceZones()[faceZones().findZoneID(staticSliderZoneName)];

        forAll (staticSliderAddr, faceI)
        {
            const face& curFace = f[staticSliderAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 0;
            }
        }
    }
}


// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::mixerFvMesh::movingPointsMask() const
{
    if (!movingPointsMaskPtr_)
    {
        calcMovingMask();
    }

    return *movingPointsMaskPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::mixerFvMesh::mixerFvMesh
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
    csPtr_
    (
        coordinateSystem::New
        (
            "coordinateSystem",
            dict_.subDict("coordinateSystem")
        )
    ),
    rpm_(readScalar(dict_.lookup("rpm"))),
    rotatingRegionMarker_
    (
        dict_.lookupOrDefault<point>("rotatingRegionMarker", csPtr_->origin())
    ),
    attachDetach_(dict_.lookupOrDefault<bool>("attachDetach", false)),
    movingPointsMaskPtr_(NULL)
{
    Info<< "Rotating region marker point: " << rotatingRegionMarker_
        << "  Attach-detach action = " << attachDetach_ << endl;

    addZonesAndModifiers();

    Info<< "Mixer mesh" << nl
        << "    origin       : " << cs().origin() << nl
        << "    axis         : " << cs().axis() << nl
        << "    rpm          : " << rpm_ << nl
        << "    marker       : " << rotatingRegionMarker_ << nl
        << "    attach-detach: " << attachDetach_
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixerFvMesh::~mixerFvMesh()
{
    deleteDemandDrivenData(movingPointsMaskPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mixerFvMesh::update()
{
    if (attached() && attachDetach_)
    {
        // Detach mesh
        Info << "Detaching rotors" << endl;
        autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

        deleteDemandDrivenData(movingPointsMaskPtr_);
    }

    // Save old points
    pointField oldPointsNew = allPoints();

    // Move points.  Rotational speed needs to be converted from rpm
    movePoints
    (
        csPtr_->globalPosition
        (
            csPtr_->localPosition(allPoints())
          + vector(0, rpm_*360.0*time().deltaT().value()/60.0, 0)
            *movingPointsMask()
        )
    );

    autoPtr<mapPolyMesh> topoChangeMap = topoChanger_.changeMesh();

    bool localMorphing = topoChangeMap->morphing();
    bool globalMorphing = localMorphing;

    reduce(globalMorphing, orOp<bool>());

    if (globalMorphing)
    {
        Info << "Attaching rotors" << endl;

        deleteDemandDrivenData(movingPointsMaskPtr_);

        // Move the sliding interface points to correct position
        if (localMorphing)
        {
            pointField mappedOldPointsNew(allPoints().size());
            mappedOldPointsNew.map(oldPointsNew, topoChangeMap->pointMap());

            movePoints(mappedOldPointsNew);

            resetMotion();
            setV0();

            // Move the sliding interface points to correct position
            movePoints(topoChangeMap->preMotionPoints());
        }
        else
        {
            pointField newPoints = allPoints();
            movePoints(oldPointsNew);

            resetMotion();
            setV0();

            // Move the sliding interface points to correct position
            movePoints(newPoints);
        }

    }

    return topoChangeMap->morphing();
}


// ************************************************************************* //
