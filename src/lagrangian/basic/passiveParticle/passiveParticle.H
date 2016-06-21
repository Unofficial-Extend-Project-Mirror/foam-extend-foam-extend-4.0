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
    Foam::passiveParticle

Description

SourceFiles
    passiveParticleI.H
    passiveParticle.C
    passiveParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef passiveParticle_H
#define passiveParticle_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class passiveParticle Declaration
\*---------------------------------------------------------------------------*/

class passiveParticle
:
    public Particle<passiveParticle>
{

public:

    // Constructors

        //- Construct from components
        passiveParticle
        (
            const Cloud<passiveParticle>& c,
            const vector& position,
            const label celli
        )
        :
            Particle<passiveParticle>(c, position, celli)
        {}

        //- Construct from Istream
        passiveParticle
        (
            const Cloud<passiveParticle>& c,
            Istream& is,
            bool readFields = true
        )
        :
            Particle<passiveParticle>(c, is, readFields)
        {}

        //- Construct as copy
        passiveParticle(const passiveParticle& p)
        :
            Particle<passiveParticle>(p)
        {}

        //- Construct and return a clone
        autoPtr<passiveParticle> clone() const
        {
            return autoPtr<passiveParticle>(new passiveParticle(*this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
