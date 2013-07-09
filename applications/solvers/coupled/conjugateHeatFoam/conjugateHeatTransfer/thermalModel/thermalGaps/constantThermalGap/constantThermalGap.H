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
    constantThermalGap

Description
    Constant thermal properties

Author
    Henrik Rusche, Wikki GmbH, All rights reserved.

SourceFiles
    constantThermalGap.C

\*---------------------------------------------------------------------------*/

#ifndef constantThermalGap_H
#define constantThermalGap_H

#include "thermalGap.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantThermalGap Declaration
\*---------------------------------------------------------------------------*/

class constantThermalGap
:
    public thermalGap
{
    // Private data

        //- thermal resistance of the gap
        dimensionedScalar beta_;

        //- list of face zones
        const wordList zones_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantThermalGap(const constantThermalGap&);

        //- Disallow default bitwise assignment
        void operator=(const constantThermalGap&);


public:

    //- Runtime type information
    TypeName("constant");

    // Static data members


    // Constructors

        //- Construct from dictionary
        constantThermalGap
        (
            const word& name,
            const volScalarField& T,
            const dictionary& dict
        );


    // Destructor

        virtual ~constantThermalGap();


    // Member Functions

        //- Return density
        virtual void modifyResistance(surfaceScalarField&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
