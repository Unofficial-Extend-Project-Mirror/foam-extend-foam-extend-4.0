/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008 H. Rusche
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

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
