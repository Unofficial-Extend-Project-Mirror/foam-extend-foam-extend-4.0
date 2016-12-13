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
    Foam::solidBodyMotionFunctions::rotatingOscillation

Description
    Rotating oscillation motion function

SourceFiles
    rotatingOscillation.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingOscillation_H
#define rotatingOscillation_H

#include "solidBodyMotionFunction.H"
#include "primitiveFields.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                      Class rotatingOscillation Declaration
\*---------------------------------------------------------------------------*/

class rotatingOscillation
:
    public solidBodyMotionFunction
{
    // Private data

        //- Centre of gravity
        point origin_;

        //- Amplitude
        vector amplitude_;

        //- Motion period
        scalar period_;


    // Private Member Functions

        //- Disallow copy construct
        rotatingOscillation(const rotatingOscillation&);

        //- Disallow default bitwise assignment
        void operator=(const rotatingOscillation&);


        //- Calculate position
        vector calcPosition(const scalar t) const;

        //- Calculate tranformation
        septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("rotatingOscillation");


    // Constructors

        //- Construct from components
        rotatingOscillation
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    //- Destructor
    virtual ~rotatingOscillation();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

    //- Return the solid-body motion velocity
        virtual septernion velocity() const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
