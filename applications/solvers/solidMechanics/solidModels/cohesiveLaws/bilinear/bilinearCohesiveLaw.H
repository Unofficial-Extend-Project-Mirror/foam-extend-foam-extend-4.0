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
    bilinearCohesiveLaw

Description
    Bilinear cohesive law.

SourceFiles
    bilinearCohesiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef bilinearCohesiveLaw_H
#define bilinearCohesiveLaw_H

#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class bilinearCohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class bilinearCohesiveLaw
:
    public cohesiveLaw
{

private:

    // Private data

        //- Transition stress
        dimensionedScalar sigma1_;

        //- Transition separation distance
        dimensionedScalar delta1_;

        //- Critical separation distance
        dimensionedScalar deltaC_;

public:

    //- Runtime type information
    TypeName("bilinear");


    // Constructors

        //- Construct from components
        bilinearCohesiveLaw
        (
            const word& cohesiveLawName,
            const dictionary& dict            
        );

        //- Construct as copy
        bilinearCohesiveLaw(const bilinearCohesiveLaw&);

        //- Construct and return a clone
        virtual autoPtr<cohesiveLaw> clone() const
        {
            return autoPtr<cohesiveLaw>(new bilinearCohesiveLaw(*this));
        }


    // Destructor

        virtual ~bilinearCohesiveLaw();


    // Member Functions

        //- Return reference to transition cohesive stress
        const dimensionedScalar& sigma1() const
        {
            return sigma1_;
        }

        //- Return reference to transition separation distance
        const dimensionedScalar& delta1() const
        {
            return delta1_;
        }

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
