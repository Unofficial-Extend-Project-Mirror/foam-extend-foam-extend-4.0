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

#include "topoBody.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChanger.H"
#include "layerAdditionRemoval.H"
#include "foamTime.H"
#include "transformField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::topoBody::addZones
(
    DynamicList<pointZone*>& pz,
    DynamicList<faceZone*>& fz,
    DynamicList<cellZone*>& cz
)
{
    // Check presence of moving cell zone
    const label mcZoneID = mesh_.cellZones().findZoneID(movingCellsName_);

    if (mcZoneID < 0)
    {
        FatalErrorIn
        (
            "void topoBody::addZones\n"
            "(\n"
            "    DynamicList<pointZone*>& pz,\n"
            "    DynamicList<faceZone*>& fz,\n"
            "    DynamicList<cellZone*>& cz\n"
            ")"
        )   << "Cannot find moving cell zone " << movingCellsName_
            << " for body " << name_
            << abort(FatalError);
    }

    forAll (layerFacesNames_, lfI)
    {
        const label fZoneID =
            mesh_.faceZones().findZoneID(layerFacesNames_[lfI]);

        if (fZoneID < 0)
        {
            FatalErrorIn
            (
                "void topoBody::addZones\n"
                "(\n"
                "    DynamicList<pointZone*>& pz,\n"
                "    DynamicList<faceZone*>& fz,\n"
                "    DynamicList<cellZone*>& cz\n"
                ")"
            )   << "Cannot find layering zone " << layerFacesNames_[lfI]
                << " for body " << name_
                << abort(FatalError);
        }
    }
}


void Foam::topoBody::addModifiers
(
    polyTopoChanger& tc,
    label& nextI
)
{
    // Add a topology modifier
    forAll (layerFacesNames_, lfI)
    {
        Info<< "Creating layering topology modifier " << layerFacesNames_[lfI]
            << " on object " << name_ << endl;

        tc.set
        (
            nextI,
            new layerAdditionRemoval
            (
                layerFacesNames_[lfI] + "Layer"  + name_,
                nextI,
                tc,
                layerFacesNames_[lfI],
                minThickness_,
                maxThickness_
            )
        );

        nextI++;
    }
}


void Foam::topoBody::calcMovingMask() const
{
    if (movingPointsMaskPtr_)
    {
        FatalErrorIn("void topoBody::calcMovingMask() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = new scalarField(mesh_.allPoints().size(), 0);
    scalarField& movingPointsMask = *movingPointsMaskPtr_;

    const cellList& c = mesh_.cells();
    const faceList& f = mesh_.allFaces();

    const labelList& cellAddr = mesh_.cellZones()
        [mesh_.cellZones().findZoneID(movingCellsName_)];

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
}


// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::topoBody::movingPointsMask() const
{
    if (!movingPointsMaskPtr_)
    {
        calcMovingMask();
    }

    return *movingPointsMaskPtr_;
}


void Foam::topoBody::clearPointMask()
{
    deleteDemandDrivenData(movingPointsMaskPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoBody::topoBody
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    movingCellsName_(dict.lookup("movingCells")),
    layerFacesNames_(dict.lookup("layerFaces")),
    minThickness_(readScalar(dict.lookup("minThickness"))),
    maxThickness_(readScalar(dict.lookup("maxThickness"))),
    SBMFPtr_(solidBodyMotionFunction::New(dict, mesh.time())),
    invertMotionMask_
    (
        dict.lookupOrDefault<bool>("invertMotionMask", false)
    ),
    movingPointsMaskPtr_(NULL)
{
    Info<< "Moving body " << name << ":" << nl
        << "    moving cells: " << movingCellsName_ << nl
        << "    layer faces : " << layerFacesNames_ << nl
        << "    invert mask : " << invertMotionMask_ << nl
        << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::topoBody::~topoBody()
{
    clearPointMask();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::topoBody::pointMotion() const
{
    // Rotational speed needs to be converted from rpm
    scalarField mpm = movingPointsMask();

    if (invertMotionMask_)
    {
        Info << "Inverting motion mask" << endl;
        mpm = 1 - mpm;
    }

    return mpm*transform(SBMFPtr_().velocity(), mesh_.allPoints())*
        mesh_.time().deltaT().value();
}


void Foam::topoBody::updateTopology()
{
    clearPointMask();
}


// ************************************************************************* //
