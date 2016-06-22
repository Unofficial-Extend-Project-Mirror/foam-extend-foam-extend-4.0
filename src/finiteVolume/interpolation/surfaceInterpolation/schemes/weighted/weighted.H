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
    Foam::weighted

Description
    Interpolation scheme class using weights looked-up from the objectRegistry.

SourceFiles
    weighted.C

\*---------------------------------------------------------------------------*/

#ifndef weighted_H
#define weighted_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class weighted Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class weighted
:
    public surfaceInterpolationScheme<Type>
{
    // Private member data

        const surfaceScalarField& weights_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const weighted&);


public:

    //- Runtime type information
    TypeName("weighted");


    // Constructors

        //- Construct from weights
        weighted(const surfaceScalarField& weights)
        :
            surfaceInterpolationScheme<Type>(weights.mesh()),
            weights_(weights)
        {}

        //- Construct from Istream
        weighted(const fvMesh& mesh, Istream& is)
        :
            surfaceInterpolationScheme<Type>(mesh),
            weights_
            (
                this->mesh().objectRegistry::template
                lookupObject<const surfaceScalarField>(word(is))
            )
        {}

        //- Construct from faceFlux and Istream
        weighted
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            weights_
            (
                this->mesh().objectRegistry::template
                lookupObject<const surfaceScalarField>
                (
                    word(is)
                )
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return weights_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
