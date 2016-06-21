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
    Foam::kinematicParcelInjectionData

Description
    Container class to provide injection data for kinematic parcels

SourceFiles
    kinematicParcelInjectionData.C

\*---------------------------------------------------------------------------*/

#ifndef kinematicParcelInjectionData_H
#define kinematicParcelInjectionData_H

#include "dictionary.H"
#include "vector.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class kinematicParcelInjectionData;

// Forward declaration of friend functions

Ostream& operator<<
(
    Ostream&,
    const kinematicParcelInjectionData&
);

Istream& operator>>
(
    Istream&,
    kinematicParcelInjectionData&
);

/*---------------------------------------------------------------------------*\
               Class kinematicParcelInjectionData Declaration
\*---------------------------------------------------------------------------*/

class kinematicParcelInjectionData
{
protected:

    // Parcel properties

        //- Position [m]
        point x_;

        //- Velocity [m/s]
        vector U_;

        //- Diameter [m]
        scalar d_;

        //- Density [kg/m3]
        scalar rho_;

        //- Mass flow rate [kg/s]
        scalar mDot_;


public:

    //- Runtime type information
    TypeName("kinematicParcelInjectionData");

    // Constructors

        //- Null constructor
        kinematicParcelInjectionData();

        //- Construct from dictionary
        kinematicParcelInjectionData(const dictionary& dict);

        //- Construct from Istream
        kinematicParcelInjectionData(Istream& is);


    //-Destructor
    virtual ~kinematicParcelInjectionData();


    // Access

        //- Return const access to the position
        inline const point& x() const;

        //- Return const access to the velocity
        inline const vector& U() const;

        //- Return const access to the diameter
        inline scalar d() const;

        //- Return const access to the density
        inline scalar rho() const;

        //- Return const access to the mass flow rate
        inline scalar mDot() const;


    // Edit

        //- Return access to the position
        inline point& x();

        //- Return access to the velocity
        inline vector& U();

        //- Return access to the diameter
        inline scalar& d();

        //- Return access to the density
        inline scalar& rho();

        //- Return access to the mass flow rate
        inline scalar& mDot();


    // I-O

        // Friend Operators

            friend bool operator==
            (
                const kinematicParcelInjectionData& a,
                const kinematicParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const kinematicParcelInjectionData&, "
                        "const kinematicParcelInjectionData&"
                    ")"
                );

                return false;
            }

            friend bool operator!=
            (
                const kinematicParcelInjectionData& a,
                const kinematicParcelInjectionData& b
            )
            {
                notImplemented
                (
                    "operator=="
                    "("
                        "const kinematicParcelInjectionData&, "
                        "const kinematicParcelInjectionData&"
                    ")"
                );

                return false;
            }

            //- Ostream operator
            friend Ostream& operator<<
            (
                Ostream& os,
                const kinematicParcelInjectionData& data
            );

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                kinematicParcelInjectionData& data
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "kinematicParcelInjectionDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
