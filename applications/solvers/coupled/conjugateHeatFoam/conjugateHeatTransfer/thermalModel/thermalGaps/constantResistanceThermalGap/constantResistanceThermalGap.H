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
