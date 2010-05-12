/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::basicReactingCloud

Description

SourceFiles
    basicReactingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef basicReactingCloud_H
#define basicReactingCloud_H

#include "ReactingCloud.H"
#include "basicReactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basicReactingCloud Declaration
\*---------------------------------------------------------------------------*/

class basicReactingCloud
:
    public ReactingCloud<basicReactingParcel>
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        basicReactingCloud(const basicReactingCloud&);

        //- Disallow default bitwise assignment
        void operator=(const basicReactingCloud&);


public:

    //- Runtime type information
    TypeName("basicReactingCloud");


    // Constructors

        //- Construct given carrier gas fields
        basicReactingCloud
        (
            const word& cloudType,
            const volPointInterpolation& vpi,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            hCombustionThermo& thermo,
            PtrList<specieReactingProperties>& gases
        );


    //- Destructor

        ~basicReactingCloud();


    // Member Functions

        //- Write fields
        void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
