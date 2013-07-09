/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 H. Rusche
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
