/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM

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

Instructions:
    This tool is used to have multiple rotating regions around the same origin
    with different rpms.
    Creating the cellZones is not implemented in this tool.
    The steps to obtain the cellZones are :

    1) use regionCellSets utility. With this command you can have different
       cellSets for each region.

    2) Change the name of the regions  eg. CellRegion0 to Rotor1 ,CellRegion1
       to Stator and vice versa.

    3) run command " setsToZones -noFlipMap ".  After this command the
       cellSets are transformed to cellZones.

    4) in dynamicMeshDict rpm section should be added.

    5) The case is ready to be run.

    Implemented by Fethi Tekin, 24.03.2010

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.

\*---------------------------------------------------------------------------*/
#include "turboFvMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(turboFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, turboFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turboFvMesh::addZonesAndModifiers()
{
    // Add zones and modifiers for motion action
    // No functionality is implemented

    if (cellZones().size() > 0)
    {
        Info<< "void turboFvMesh::addZonesAndModifiers() : "
            << "Zones and modifiers already present.  Skipping."
            << endl;

        return;
    }
    else
    {
            FatalErrorIn("turboFvMesh")
            << "Cell Regions have to be created"
            << abort(FatalError);
    }
}
void Foam::turboFvMesh::calcMovingPoints() const
{
    if (debug)
    {
        Info<< "void turboFvMesh::calcMovingMasks() const : "
            << "Calculating point and cell masks"
            << endl;
    }

    if (movingPointsPtr_)
    {
        FatalErrorIn("void turboFvMesh::calcMovingMasks() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Retrieve the zone Names
    const wordList zoneNames = cellZones().names();

    // Set the points
    movingPointsPtr_ = new vectorField(allPoints().size(),vector::zero);

    vectorField& movingPoints = *movingPointsPtr_;

    const cellList& c = cells();
    const faceList& f = allFaces();

    forAll(zoneNames,zoneI)
    {
        Info<< "Moving Region Zone Name:" << zoneNames[zoneI]<< nl<<endl;

        const labelList& cellAddr =
        cellZones()[cellZones().findZoneID(zoneNames[zoneI])];

        rpm_ = readScalar(dict_.subDict("rpm").lookup(zoneNames[zoneI]));

        Info<< "rpm:" << rpm_<<nl<<endl;

        forAll (cellAddr, cellI)
        {
            const cell& curCell = c[cellAddr[cellI]];

            forAll (curCell, faceI)
            {
                 const face& curFace = f[curCell[faceI]];

                 forAll (curFace, pointI)
                 {
                     //The rotation data is saved within the cell data. For
                     //non-rotating regions rpm is zero,so mesh movement is
                     //also zero. The conversion of rotational speed
                     movingPoints[curFace[pointI]] =
                         vector(0,rpm_*360.0*time().deltaT().value()/60.0, 0);
                 }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::turboFvMesh::turboFvMesh
(
    const IOobject& io
)
:
    dynamicFvMesh(io),
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

    movingPointsPtr_(NULL)
{
    addZonesAndModifiers();

    Info<< "Turbomachine Mixer mesh:" << nl
        << "    origin: " << cs().origin() << nl
        << "    axis  : " << cs().axis() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turboFvMesh::~turboFvMesh()
{
    deleteDemandDrivenData(movingPointsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return moving points mask
const Foam::vectorField& Foam::turboFvMesh::movingPoints() const
{
    if (!movingPointsPtr_)
    {
        calcMovingPoints();
    }

    return *movingPointsPtr_;
}


bool Foam::turboFvMesh::update()
{
    movePoints
    (
        csPtr_->globalPosition
        (
          csPtr_->localPosition(allPoints()) + movingPoints()
        )
    );

    // The mesh is not morphing
    return false;
}


// ************************************************************************* //
