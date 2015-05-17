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

#include "constantThermalSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantThermalSource, 0);
    addToRunTimeSelectionTable(thermalSource, constantThermalSource, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::constantThermalSource::constantThermalSource
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict
)
:
    thermalSource(name, T, dict),
    S_(dict.lookup("S")),
    zones_(dict.lookup("zones"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantThermalSource::~constantThermalSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constantThermalSource::addSource(volScalarField& source) const
{
    forAll(zones_, zoneI)
    {
        const label zoneID = mesh().cellZones().findZoneID(zones_[zoneI]);

        if ( zoneID < 0 )
        {
            FatalErrorIn
            (
                "constantThermalSource::addSourcex()\n"
            )   << "Zone " << zones_[zoneI]
                << " specified in source " << name()
                << " does not exist"
                << abort(FatalError);
        }

        const labelList& cells = mesh().cellZones()[zoneID];

        forAll(cells, cellI)
        {
            source[cells[cellI]] = S_.value();
        }
    }
}


// ************************************************************************* //
