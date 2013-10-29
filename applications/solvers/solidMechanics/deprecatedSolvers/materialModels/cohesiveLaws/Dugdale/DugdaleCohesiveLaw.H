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
    DugdaleCohesiveLaw

Description
    Dugdale cohesive law.

SourceFiles
    DugdaleCohesiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef DugdaleCohesiveLaw_H
#define DugdaleCohesiveLaw_H

#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class DugdaleCohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class DugdaleCohesiveLaw
:
    public cohesiveLaw
{
private:

    // Private data

        //- Critical separation distance
        dimensionedScalar deltaC_;


public:

    //- Runtime type information
    TypeName("Dugdale");


    // Constructors

        //- Construct from components
        DugdaleCohesiveLaw
        (
            const word& cohesiveLawName,
            const dictionary& dict            
        );

        //- Construct as copy
        DugdaleCohesiveLaw(const DugdaleCohesiveLaw&);

        //- Construct and return a clone
        virtual autoPtr<cohesiveLaw> clone() const
        {
            return autoPtr<cohesiveLaw>(new DugdaleCohesiveLaw(*this));
        }


    // Destructor

        virtual ~DugdaleCohesiveLaw();


    // Member Functions

        //- Return reference to critical separation distance
        const dimensionedScalar& deltaC() const
        {
            return deltaC_;
        }

        //- Return current holding traction
        scalar traction(scalar delta) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
