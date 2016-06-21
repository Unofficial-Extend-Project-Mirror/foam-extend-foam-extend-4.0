/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    pressureVaryingCoulomb

Description
    pressureVaryingCoulomb friction law.

    Shear traction is linearly related to the contact pressure via the friction
    coefficient, where the friction coefficient is specified as a time series.

SourceFiles
    pressureVaryingCoulomb.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef pressureVaryingCoulomb_H
#define pressureVaryingCoulomb_H

#include "frictionLaw.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class pressureVaryingCoulomb Declaration
\*---------------------------------------------------------------------------*/

class pressureVaryingCoulomb
:
    public frictionLaw
{
    // Private data

         //- friction law coefficientdict
        dictionary frictionLawDict_;

        //- Coulomb friction coefficient time series
        interpolationTable<scalar> frictionCoeffSeries_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        //pressureVaryingCoulomb(const pressureVaryingCoulomb&);

        //- Disallow default bitwise assignment
        void operator=(const pressureVaryingCoulomb&);

public:

    //- Runtime type information
    TypeName("pressureVaryingCoulomb");

    // Static data members


    // Constructors

        //- Construct from dictionary
        pressureVaryingCoulomb
        (
            const word& name,
            const frictionContactModel& fricModel,
            const dictionary& dict
        );

        //- Construcct as a copy
        pressureVaryingCoulomb(const pressureVaryingCoulomb& fricLaw);

        //- Construct and return a clone
        virtual autoPtr<frictionLaw> clone() const
        {
            return
                autoPtr<frictionLaw>
                (
                    new pressureVaryingCoulomb(*this)
                );
        }



    // Destructor

        virtual ~pressureVaryingCoulomb();


    // Member Functions

        //- Return slip traction field
        // the traction that would cause slipping
        // for the given pressure
        virtual scalar slipTraction(const scalar pressure);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
