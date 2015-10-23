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
    constantResistanceThermalGap

Description
    Constant thermal properties

Author
    Henrik Rusche, Wikki GmbH, All rights reserved.

SourceFiles
    constantResistanceThermalGap.C

\*---------------------------------------------------------------------------*/

#ifndef constantResistanceThermalGap_H
#define constantResistanceThermalGap_H

#include "thermalGap.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantResistanceThermalGap Declaration
\*---------------------------------------------------------------------------*/

class constantResistanceThermalGap
:
    public thermalGap
{
    // Private data

        //- thermal resistance of the gap
        dimensionedScalar R_;

        //- list of face zones
        const wordList zones_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantResistanceThermalGap(const constantResistanceThermalGap&);

        //- Disallow default bitwise assignment
        void operator=(const constantResistanceThermalGap&);


public:

    //- Runtime type information
    TypeName("constantResistance");

    // Static data members


    // Constructors

        //- Construct from dictionary
        constantResistanceThermalGap
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~constantResistanceThermalGap();


    // Member Functions

        //- Return density
        virtual void modifyResistance(surfaceScalarField&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
