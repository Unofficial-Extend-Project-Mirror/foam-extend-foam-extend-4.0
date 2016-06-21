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
    Foam::reverseLinear

Description
    Inversed weight central-differencing interpolation scheme class.

    Useful for inverse weighted and harmonic interpolations.

SourceFiles
    reverseLinear.C

\*---------------------------------------------------------------------------*/

#ifndef reverseLinear_H
#define reverseLinear_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reverseLinear Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class reverseLinear
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const reverseLinear&);


public:

    //- Runtime type information
    TypeName("reverseLinear");


    // Constructors

        //- Construct from mesh
        reverseLinear(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        reverseLinear(const fvMesh& mesh, Istream&)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        reverseLinear
        (
            const fvMesh& mesh,
            const surfaceScalarField&,
            Istream&
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            const fvMesh& mesh = this->mesh();

            tmp<surfaceScalarField> tcdWeights
            (
                mesh.surfaceInterpolation::weights()
            );
            const surfaceScalarField& cdWeights = tcdWeights();

            tmp<surfaceScalarField> treverseLinearWeights
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "reverseLinearWeights",
                        mesh.time().timeName(),
                        mesh
                    ),
                    mesh,
                    dimless
                )
            );
            surfaceScalarField& reverseLinearWeights = treverseLinearWeights();

            reverseLinearWeights.internalField() =
                1.0 - cdWeights.internalField();

            forAll (mesh.boundary(), patchI)
            {
                if (mesh.boundary()[patchI].coupled())
                {
                    reverseLinearWeights.boundaryField()[patchI] =
                        1.0 - cdWeights.boundaryField()[patchI];
                }
                else
                {
                    reverseLinearWeights.boundaryField()[patchI] =
                        cdWeights.boundaryField()[patchI];
                }
            }

            return treverseLinearWeights;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
