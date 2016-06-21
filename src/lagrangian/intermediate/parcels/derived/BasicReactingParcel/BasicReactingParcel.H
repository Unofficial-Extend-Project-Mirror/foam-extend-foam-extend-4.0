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
    Foam::BasicReactingParcel

Description


SourceFiles
    BasicReactingParcel.C
    BasicReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BasicReactingParcel_H
#define BasicReactingParcel_H

#include "objectRegistry.H"
#include "ReactingParcelTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class ThermoType>
class BasicReactingParcel;

/*---------------------------------------------------------------------------*\
                    Class BasicReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class BasicReactingParcel
:
    public ReactingParcel<BasicReactingParcel<ThermoType> >
{

public:

    //- The type of thermodynamics this parcel was instantiated for
    typedef ThermoType thermoType;

    //- Run-time type information
    TypeName("BasicReactingParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        BasicReactingParcel
        (
            ReactingCloud<BasicReactingParcel>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        BasicReactingParcel
        (
            ReactingCloud<BasicReactingParcel>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const scalarField& Y0,
            const typename ReactingParcel<BasicReactingParcel>::
                constantProperties& constProps
        );

        //- Construct from Istream
        BasicReactingParcel
        (
            const Cloud<BasicReactingParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        BasicReactingParcel(const BasicReactingParcel& p);

        //- Construct and return a clone
        autoPtr<BasicReactingParcel> clone() const
        {
            return
                autoPtr<BasicReactingParcel>
                (
                    new BasicReactingParcel(*this)
                );
        }


    //- Destructor
    virtual ~BasicReactingParcel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
