/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "constantResistanceThermalGap.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantResistanceThermalGap, 0);
    addToRunTimeSelectionTable(thermalGap, constantResistanceThermalGap, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::constantResistanceThermalGap::constantResistanceThermalGap
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict
)
:
    thermalGap(name, T, dict),
    R_(dict.lookup("R")),
    zones_(dict.lookup("zones"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantResistanceThermalGap::~constantResistanceThermalGap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantResistanceThermalGap::modifyResistance
(
    surfaceScalarField& DT
) const
{
    forAll(zones_, zoneI)
    {
        const label zoneID = mesh().faceZones().findZoneID(zones_[zoneI]);

        if ( zoneID < 0 )
        {
            FatalErrorIn
            (
                "constantResistanceThermalGap::modifyResistance()\n"
            )   << "Zone " << zones_[zoneI]
                << " specified in gap " << name()
                << " does not exist"
                << abort(FatalError);
        }

        const labelList& faces = mesh().faceZones()[zoneID];
        const surfaceScalarField& deltaCoeffs = mesh().deltaCoeffs();

        forAll(faces, faceI)
        {
            DT[faces[faceI]] = 1.0/deltaCoeffs[faceI]/R_.value();
        }
    }
}


// ************************************************************************* //
