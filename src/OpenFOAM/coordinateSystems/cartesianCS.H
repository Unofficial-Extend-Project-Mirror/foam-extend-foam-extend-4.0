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
    Foam::cartesianCS

Description
    General coordinate transformations from arbitrary coordinate systems
    to the global Cartesian system

SourceFiles
    cartesianCS.C

\*---------------------------------------------------------------------------*/

#ifndef cartesianCS_H
#define cartesianCS_H

#include "coordinateSystem.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cartesianCS Declaration
\*---------------------------------------------------------------------------*/

class cartesianCS
:
    public coordinateSystem
{
public:

    //- Runtime type information
    TypeName("cartesian");


    // Constructors

        //- Construct null
        cartesianCS();

        //- Construct from origin and 2 axes
        cartesianCS
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& dir
        );

        //- Construct from origin and rotation
        cartesianCS
        (
            const word& name,
            const point& origin,
            const coordinateRotation& cr
        );

        //- Construct from dictionary
        cartesianCS(const word& name, const dictionary& dict);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
