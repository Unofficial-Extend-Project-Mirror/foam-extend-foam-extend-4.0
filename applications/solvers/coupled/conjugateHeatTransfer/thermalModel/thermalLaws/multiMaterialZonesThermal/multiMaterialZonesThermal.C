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

#include "multiMaterialZonesThermal.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiMaterialZonesThermal, 0);
    addToRunTimeSelectionTable(thermalLaw, multiMaterialZonesThermal, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiMaterialZonesThermal::multiMaterialZonesThermal
(
    const word& name,
    const volScalarField& T,
    const dictionary& dict
)
:
    multiMaterialThermal(name, T, dict, -1)
{
    readLaws(T, dict);

    PtrList<entry> lawEntries(dict.lookup("laws"));

    PtrList<thermalLaw>& laws = *this;
    forAll (laws, lawI)
    {
        wordList zones (lawEntries[lawI].dict().lookup("zones"));

        forAll(zones, zoneI)
        {
            const label zoneID = mesh().cellZones().findZoneID(zones[zoneI]);

            if ( zoneID < 0 )
            {
                FatalErrorIn
                (
                    "multiMaterialZonesThermal::multiMaterialZonesThermal()\n"
                )   << "Zone " << zones[zoneI]
                    << " specified in material " << lawEntries[lawI].keyword()
                    << " does not exist"
                    << abort(FatalError);
            }

            const labelList& cells = mesh().cellZones()[zoneID];

            forAll(cells, cellI)
            {
                materials_[cells[cellI]] = lawI;
            }
        }
    }

    materials_.correctBoundaryConditions();
    materials_.write();

    checkLaws();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiMaterialZonesThermal::~multiMaterialZonesThermal()
{}

// ************************************************************************* //
