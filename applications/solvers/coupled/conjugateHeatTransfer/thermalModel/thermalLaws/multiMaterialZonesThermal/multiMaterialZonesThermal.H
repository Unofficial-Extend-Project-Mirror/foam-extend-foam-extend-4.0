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
