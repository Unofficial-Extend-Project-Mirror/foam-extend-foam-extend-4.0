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

#include "mixerRotor.H"
#include "regionSplit.H"
#include "polyTopoChanger.H"
#include "slidingInterface.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mixerRotor::addZones
(
    DynamicList<pointZone*>& pz,
    DynamicList<faceZone*>& fz,
    DynamicList<cellZone*>& cz,
    const regionSplit& rs
)
{
    // Get the region of the cell containing the origin.
    label originRegion = rs[mesh_.findNearestCell(rotatingRegionMarker_)];

    labelList movingCells(mesh_.nCells());
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
    Info<< "Number of moving cells for " << name_ << ": "
        << nMovingCells << endl;

    cz.append
    (
        new cellZone
        (
            "movingCellsZone" + name_,
            movingCells,
            cz.size(),
            mesh_.cellZones()
        )
    );

    if (useTopoSliding_)
    {
        // Add an empty zone for cut points
        pz.append
        (
            new pointZone
            (
                "cutPointZone" + name_,
                labelList(0),
                pz.size(),
                mesh_.pointZones()
            )
        );


        // Do face zones for slider

        // Moving slider
        label movingSliderIndex =
            mesh_.boundaryMesh().findPatchID(movingSliderName_);

        if (movingSliderIndex < 0)
        {
            FatalErrorIn("void mixerRotor::addZones(...) const")
                << "Moving slider patch not found in boundary"
                << abort(FatalError);
        }

        label staticSliderIndex =
            mesh_.boundaryMesh().findPatchID(staticSliderName_);

        if (staticSliderIndex < 0)
        {
            FatalErrorIn("void mixerRotor::addZones(...) const")
                << "Static slider patch not found in boundary"
                << abort(FatalError);
        }

        const polyPatch& movingSlider =
            mesh_.boundaryMesh()[movingSliderIndex];

        labelList isf(movingSlider.size());

        forAll (isf, i)
        {
            isf[i] = movingSlider.start() + i;
        }

        fz.append
        (
            new faceZone
            (
                movingSliderName_ + "Zone" + name_,
                isf,
                boolList(movingSlider.size(), false),
                fz.size(),
                mesh_.faceZones()
            )
        );

        // Static slider
        const polyPatch& staticSlider =
            mesh_.boundaryMesh()[staticSliderIndex];

        labelList osf(staticSlider.size());

        forAll (osf, i)
        {
            osf[i] = staticSlider.start() + i;
        }

        fz.append
        (
            new faceZone
            (
                staticSliderName_ + "Zone" + name_,
                osf,
                boolList(staticSlider.size(), false),
                fz.size(),
                mesh_.faceZones()
            )
        );

        // Add empty zone for cut faces
        fz.append
        (
            new faceZone
            (
                "cutFaceZone" + name_,
                labelList(0),
                boolList(0, false),
                fz.size(),
                mesh_.faceZones()
            )
        );
    }
}


void Foam::mixerRotor::addModifiers
(
    polyTopoChanger& tc,
    label& nextI
)
{
    // Add a topology modifier
    if (useTopoSliding())
    {
        Info << "Adding topology modifier for rotor " << name_ << endl;

        tc.set
        (
            nextI,
            new slidingInterface
            (
                "mixerSlider"  + name_,
                nextI,
                tc,
                staticSliderName_ + "Zone"  + name_,
                movingSliderName_ + "Zone" + name_,
                "cutPointZone" + name_,
                "cutFaceZone" + name_,
                staticSliderName_,
                movingSliderName_,
                slidingInterface::INTEGRAL,   // Edge matching algorithm
                attachDetach_,                // Attach-detach action
                intersection::VISIBLE         // Projection algorithm
            )
        );

        nextI++;
    }
}


void Foam::mixerRotor::calcMovingMask() const
{
    if (movingPointsMaskPtr_)
    {
        FatalErrorIn("void mixerRotor::calcMovingMask() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = new scalarField(mesh_.allPoints().size(), 0);
    scalarField& movingPointsMask = *movingPointsMaskPtr_;

    const cellList& c = mesh_.cells();
    const faceList& f = mesh_.allFaces();

    const labelList& cellAddr = mesh_.cellZones()
        [mesh_.cellZones().findZoneID("movingCellsZone" + name_)];

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

    // Attempt to enforce motion on sliders if zones exist
    const label msI =
        mesh_.faceZones().findZoneID(movingSliderName_ + "Zone" + name_);

    if (msI > -1)
    {
        const labelList& movingSliderAddr = mesh_.faceZones()[msI];

        forAll (movingSliderAddr, faceI)
        {
            const face& curFace = f[movingSliderAddr[faceI]];

            forAll (curFace, pointI)
            {
                movingPointsMask[curFace[pointI]] = 1;
            }
        }
    }

    const label ssI =
        mesh_.faceZones().findZoneID(staticSliderName_ + "Zone" + name_);

    if (ssI > -1)
    {
        const labelList& staticSliderAddr = mesh_.faceZones()[ssI];

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
const Foam::scalarField& Foam::mixerRotor::movingPointsMask() const
{
    if (!movingPointsMaskPtr_)
    {
        calcMovingMask();
    }

    return *movingPointsMaskPtr_;
}


void Foam::mixerRotor::clearPointMask()
{
    deleteDemandDrivenData(movingPointsMaskPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixerRotor::mixerRotor
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    csPtr_
    (
        coordinateSystem::New
        (
            "coordinateSystem",
            dict.subDict("coordinateSystem")
        )
    ),
    rpm_(readScalar(dict.lookup("rpm"))),
    movingSliderName_(dict.lookup("movingPatch")),
    staticSliderName_(dict.lookup("staticPatch")),
    rotatingRegionMarker_
    (
        dict.lookupOrDefault<point>("rotatingRegionMarker", csPtr_->origin())
    ),
    invertMotionMask_
    (
        dict.lookupOrDefault<bool>("invertMotionMask", false)
    ),
    useTopoSliding_(dict.lookup("useTopoSliding")),
    attachDetach_(dict.lookupOrDefault<bool>("attachDetach", true)),
    movingPointsMaskPtr_(NULL)
{
    Info<< "Rotor " << name << ":" << nl
        << "    origin      : " << cs().origin() << nl
        << "    axis        : " << cs().axis() << nl
        << "    rpm         : " << rpm_ << nl
        << "    topo sliding: " << useTopoSliding_ << endl;

    if (useTopoSliding_)
    {
        Info<< "    attach-detach: " << attachDetach_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixerRotor::~mixerRotor()
{
    clearPointMask();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::mixerRotor::pointMotion() const
{
    // Rotational speed needs to be converted from rpm
    scalarField mpm = movingPointsMask();

    if (invertMotionMask_)
    {
        Info << "Inverting motion mask" << endl;
        mpm = 1 - mpm;
    }

    return csPtr_->globalPosition
    (
        // Motion vector in cymindrical coordinate system (x theta z)
        csPtr_->localPosition(mesh_.allPoints())
      + vector(0, rpm_*360.0*mesh_.time().deltaT().value()/60.0, 0)*mpm
    ) - mesh_.allPoints();
}


void Foam::mixerRotor::updateTopology()
{
    clearPointMask();
}


// ************************************************************************* //
