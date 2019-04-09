/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::decayingVorton

Description
    Foam::decayingVorton

SourceFiles
    decayingVorton.C

Author
    Prof. Nikolai Kornev, Uni Rostock
    Clean-up by Hrvoje Jasak

\*---------------------------------------------------------------------------*/

#ifndef decayingVorton_H
#define decayingVorton_H

#include "vector.H"
#include "IOstreams.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class decayingVorton;

bool operator==(const decayingVorton&, const decayingVorton&);
bool operator!=(const decayingVorton&, const decayingVorton&);

Ostream& operator<<(Ostream& s, const decayingVorton& vt);
Istream& operator>>(Istream& s, decayingVorton& vt);


/*---------------------------------------------------------------------------*\
                        Class decayingVorton Declaration
\*---------------------------------------------------------------------------*/

class decayingVorton
{
    // Private data

        //- Length-scale for a vorton
        scalar length_;

        //- Location of vorton centroid
        vector location_;

        //- Translational velocity
        vector velocity_;

        //- Max distance
        scalar xmax_;

        //- Rotational velocity
        vector omega_;


    // Static Private Members

        //- Random number generator for vortons
        static Random ranGen;


public:

    // Constructors

        //- Construct null
        decayingVorton()
        {}

        //- Construct from components
        decayingVorton
        (
            const scalar length,
            const vector location,
            const vector velocity,
            const scalar xmax
        );

        //- Construct from Istream
        decayingVorton(Istream& is);


    // Destructor - default


    // Member Functions

        // Access

            //- Return location
            const vector& location() const
            {
                return location_;
            }

            //- Return velocity
            const vector& velocity() const
            {
                return velocity_;
            }

            //- Return length
            scalar length() const
            {
                return length_;
            }

            //- Return lspot
            scalar lspot() const
            {
                return 3*length_;
            }

            //- Return xmax
            scalar xmax() const
            {
                return xmax_;
            }

            //- Return velocity at given point
            inline vector velocityAt(const vector& pos) const
            {
                vector dv = pos - location_;
                scalar nrm2 = magSqr(dv);
                vector t = dv ^ omega_;

                return (1/length_)*
                    exp
                    (
                        -(mathematicalConstant::pi/2)*nrm2/
                        (length_*length_)
                    )*t;
            }


        // Edit

            //- Move vorton
            void move(scalar dt)
            {
                location_ += vector(velocity_.x(), 0, 0)*dt;
            }


    // Friend Operators

        friend bool operator==(const decayingVorton&, const decayingVorton&);

        friend bool operator!=(const decayingVorton&, const decayingVorton&);


    // IOstream Operators

        friend Ostream& operator<<(Ostream& s, const decayingVorton& vt);

        friend Istream& operator>>(Istream& s, decayingVorton& vt);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
