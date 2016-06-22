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
    Foam::solidBodyMotionFunctions::SDA

Description
    Ship design analysis (SDA) 3DoF motion function.

    Comprising sinusoidal roll (rotation about x), heave (z-translation)
    and sway (y-translation) motions with changing amplitude and phase.

See Also
    SKA (Sea Keeping Analysis) for 6DoF motion.

SourceFiles
    SDA.C

\*---------------------------------------------------------------------------*/

#ifndef SDA_H
#define SDA_H

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class SDA Declaration
\*---------------------------------------------------------------------------*/

class SDA
:
    public solidBodyMotionFunction
{
    // Private data

        //- Center of gravity
        vector CofG_;

        //- Model scale ratio
        scalar lamda_;

        //- Max roll amplitude [rad]
        scalar rollAmax_;

        //- Min roll amplitude [rad]
        scalar rollAmin_;

        //- Heave amplitude [m]
        scalar heaveA_;

        //- Sway amplitude [m]
        scalar swayA_;

        //- Damping Coefficient [-]
        scalar Q_;

        //- Time Period for liquid [sec]
        scalar Tp_;

        //- Natural Period of Ship [sec]
        scalar Tpn_;

        //- Reference time step [sec]
        scalar dTi_;

        //- Incr. in Tp/unit 'dTi'[-]
        scalar dTp_;


    // Private Member Functions

        //- Disallow copy construct
        SDA(const SDA&);

        //- Disallow default bitwise assignment
        void operator=(const SDA&);


        //- Calculate transform given time
        septernion calcTransformation(const scalar t) const;


public:

    //- Runtime type information
    TypeName("SDA");


    // Constructors

        //- Construct from components
        SDA
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );


    // Destructor

        virtual ~SDA();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation() const;

        //- Return the solid-body motion velocity septernion
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
