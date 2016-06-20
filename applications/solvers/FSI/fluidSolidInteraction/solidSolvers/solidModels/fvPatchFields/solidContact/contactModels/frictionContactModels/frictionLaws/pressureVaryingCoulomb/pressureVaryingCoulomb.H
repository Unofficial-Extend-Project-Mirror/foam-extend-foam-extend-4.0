/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
