/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "ggiRotor.H"
#include "foamTime.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ggiRotor::calcMovingMask() const
{
    if (movingPointsMaskPtr_)
    {
        FatalErrorIn("void ggiRotor::calcMovingMask() const")
            << "point mask already calculated"
            << abort(FatalError);
    }

    // Set the point mask
    movingPointsMaskPtr_ = new scalarField(mesh_.allPoints().size(), 0);
    scalarField& movingPointsMask = *movingPointsMaskPtr_;

    const cellList& c = mesh_.cells();
    const faceList& f = mesh_.allFaces();

    const labelList& cellAddr = mesh_.cellZones()
        [mesh_.cellZones().findZoneID(movingCellsZoneName_)];

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
const Foam::scalarField& Foam::ggiRotor::movingPointsMask() const
{
    if (!movingPointsMaskPtr_)
    {
        calcMovingMask();
    }

    return *movingPointsMaskPtr_;
}


void Foam::ggiRotor::clearPointMask()
{
    deleteDemandDrivenData(movingPointsMaskPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiRotor::ggiRotor
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    cs_
    (
        "coordinateSystem",
        dict.subDict("coordinateSystem")
    ),
    rpm_(readScalar(dict.lookup("rpm"))),
    movingCellsZoneName_(dict.lookup("movingCells")),
    movingPointsMaskPtr_(NULL)
{
    // Make sure the coordinate system does not operate in degrees
    // Bug fix, HJ, 3/Oct/2011
    if (!cs_.inDegrees())
    {
        WarningIn
        (
            "ggiRotor::ggiRotor\n"
            "(\n"
            "    const word& name,\n"
            "    const polyMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")"
        )   << "Mixer coordinate system is set to operate in radians.  "
            << "Changing to rad for correct calculation of angular velocity."
            << nl
            << "To remove this message please add entry" << nl << nl
            << "inDegrees true;" << nl << nl
            << "to the specification of the coordinate system"
            << endl;

        cs_.inDegrees() = true;
    }

    Info<< "Rotor " << name << ":" << nl
        << "    origin      : " << cs().origin() << nl
        << "    axis        : " << cs().axis() << nl
        << "    rpm         : " << rpm_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ggiRotor::~ggiRotor()
{
    clearPointMask();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::ggiRotor::pointMotion() const
{
    // Rotational speed needs to be converted from rpm
    scalarField mpm = movingPointsMask();

    return cs_.globalPosition
    (
        // Motion vector in cylindrical coordinate system (x theta z)
        cs_.localPosition(mesh_.allPoints())
      + vector(0, rpm_*360.0*mesh_.time().deltaT().value()/60.0, 0)*mpm
    ) - mesh_.allPoints();
}


void Foam::ggiRotor::updateTopology()
{
    clearPointMask();
}


// ************************************************************************* //
