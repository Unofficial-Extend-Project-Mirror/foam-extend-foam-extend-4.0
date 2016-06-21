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
    Foam::TimeState

Description
    The time value with time-stepping information, user-defined remapping, etc.

SourceFiles
    TimeState.C

\*---------------------------------------------------------------------------*/

#ifndef TimeState_H
#define TimeState_H

#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TimeState Declaration
\*---------------------------------------------------------------------------*/


class TimeState
:
    public dimensionedScalar
{

protected:

        label timeIndex_;
        scalar deltaT_;
        scalar deltaTSave_;
        scalar deltaT0_;
        bool deltaTchanged_;

        label outputTimeIndex_;
        bool  outputTime_;


public:

    // Constructors

        TimeState();


    // Destructor

        //- Virtual destructor
        virtual ~TimeState();


    // Member functions

        // Access

            //- Convert the user-time (e.g. CA deg) to real-time (s).
            virtual scalar userTimeToTime(const scalar theta) const;

            //- Convert the real-time (s) into user-time (e.g. CA deg)
            virtual scalar timeToUserTime(const scalar t) const;

            //- Return current time value
            virtual scalar timeOutputValue() const;

            //- Return current time index
            virtual label timeIndex() const;

            //- Return time step value
            inline scalar deltaTValue() const
            {
                return deltaT_;
            }

            //- Return old time step value
            inline scalar deltaT0Value() const
            {
                return deltaT0_;
            }

            //- Return time step
            virtual dimensionedScalar deltaT() const;

            //- Return old time step
            virtual dimensionedScalar deltaT0() const;


        // Check

            //- Return true if this is an output time
            virtual bool outputTime() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
