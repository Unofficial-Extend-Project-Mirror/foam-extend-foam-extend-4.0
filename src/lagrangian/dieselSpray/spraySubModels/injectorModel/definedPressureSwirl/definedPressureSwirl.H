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
    Foam::definedPressureSwirlInjector

Description
    Pressure swirl injector model with user-defined
    cone-angle/interval vs. time

\*---------------------------------------------------------------------------*/

#ifndef definedPressureSwirlInjector_H
#define definedPressureSwirlInjector_H

#include "injectorModel.H"
#include "scalarList.H"
#include "vector.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class definedPressureSwirlInjector Declaration
\*---------------------------------------------------------------------------*/

class definedPressureSwirlInjector
:
    public injectorModel
{

private:

    // Private data

        typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

        dictionary definedPressureSwirlInjectorDict_;

        List<pair> coneAngle_;
        List<pair> coneInterval_;

        List<pair>  maxKv_;

    // The initial spray angle for the parcels
        mutable scalar angle_;

    // The initial velocity for the parcels
        mutable scalar u_;

    // private member functions

        scalar kv
        (
            const label inj,
            const scalar massFlow,
            const scalar dPressure,
            const scalar t
        ) const;

        scalar deltaPressureInj(const scalar time, const label inj) const;
        scalar averagePressure(const label inj) const;

public:

    //- Runtime type information
        TypeName("definedPressureSwirlInjector");


    // Constructors

        //- Construct from components
        definedPressureSwirlInjector
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~definedPressureSwirlInjector();


    // Member Functions

        //- Return the injected droplet diameter
        scalar d0(const label injector, const scalar time) const;

        //- Return the spray angle of the injector
        vector direction
        (
            const label injector,
            const label hole,
            const scalar time,
            const scalar d
        ) const;

        scalar velocity
        (
            const label i,
            const scalar time
        ) const;

        scalar averageVelocity
        (
            const label i
        ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
