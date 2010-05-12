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
    Foam::sphericalCS

Description
    Spherical coordinate system

SourceFiles
    sphericalCS.C

\*---------------------------------------------------------------------------*/

#ifndef sphericalCS_H
#define sphericalCS_H

#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class sphericalCS Declaration
\*---------------------------------------------------------------------------*/

class sphericalCS
:
    public coordinateSystem
{

protected:

    // Protected Member Functions

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual vector localToGlobal(const vector&, bool translate) const;

        //- Convert from local coordinate system to the global Cartesian system
        //  with optional translation for the origin
        virtual tmp<vectorField> localToGlobal
        (
            const vectorField&,
            bool translate
        ) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual vector globalToLocal(const vector&, bool translate) const;

        //- Convert from global Cartesian system to the local coordinate system
        //  with optional translation for the origin
        virtual tmp<vectorField> globalToLocal
        (
            const vectorField&,
            bool translate
        ) const;

public:

    //- Runtime type information
    TypeName("spherical");


    // Constructors

        //- Construct null
        sphericalCS();

        //- Construct from components
        sphericalCS
        (
            const word& name,
            const point& origin,
            const vector& axis,
            const vector& direction
        );

        //- Construct from origin and rotation angles
        sphericalCS
        (
            const word& name,
            const point& origin,
            const coordinateRotation& cr
        );


        //- Construct from dictionary
        sphericalCS(const word& name, const dictionary& dict);


    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
