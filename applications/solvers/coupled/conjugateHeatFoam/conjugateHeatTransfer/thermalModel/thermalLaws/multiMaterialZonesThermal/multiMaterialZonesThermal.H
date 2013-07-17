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

Class
    multiMaterialZonesThermal

Description
    Zoned multi-material thermal model constructed from zones

Author
    Henrik Rusche, Wikki GmbH, All rights reserved.

SourceFiles
    multiMaterialZonesThermal.C

\*---------------------------------------------------------------------------*/

#ifndef multiMaterialZonesThermal_H
#define multiMaterialZonesThermal_H

#include "multiMaterialThermal.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multiMaterialZonesThermal Declaration
\*---------------------------------------------------------------------------*/

class multiMaterialZonesThermal
:
    public multiMaterialThermal
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiMaterialZonesThermal(const multiMaterialZonesThermal&);

        //- Disallow default bitwise assignment
        void operator=(const multiMaterialZonesThermal&);

public:

    //- Runtime type information
    TypeName("multiMaterialZones");

    // Static data members


    // Constructors

        //- Construct from dictionary
        multiMaterialZonesThermal
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~multiMaterialZonesThermal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
