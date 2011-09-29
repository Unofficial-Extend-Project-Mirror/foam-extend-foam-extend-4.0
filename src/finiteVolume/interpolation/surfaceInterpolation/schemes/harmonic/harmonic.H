/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::harmonic

Description
    Harmonic-mean differencing scheme class.

    This scheme interpolates 1/field using a scheme specified at run-time
    and return the reciprocal of the interpolate.

SourceFiles
    harmonic.C

\*---------------------------------------------------------------------------*/

#ifndef harmonic_H
#define harmonic_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "magLongDelta.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class harmonic Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class harmonic
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const harmonic&);


public:

    //- Runtime type information
    TypeName("harmonic");


    // Constructors

        //- Construct from mesh
        harmonic(const fvMesh& mesh)
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from Istream
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        harmonic
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}

        //- Construct from faceFlux and Istream
        harmonic
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh)
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi
        ) const
        {
            tmp<surfaceScalarField> tw
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "harmonicWeightingFactors" + phi.name(),
                        this->mesh().time().timeName(),
                        this->mesh()
                    ),
                    this->mesh() ,
                    dimless
                )
            );

            const surfaceScalarField& deltaCoeffs = this->mesh().deltaCoeffs();
            const surfaceScalarField& weights = this->mesh().weights();

            const magLongDelta& mld = magLongDelta::New(this->mesh());
            const surfaceScalarField& longDelta = mld.magDelta();

            surfaceScalarField& w = tw();

            const unallocLabelList& owner = this->mesh().owner();
            const unallocLabelList& neighbour = this->mesh().neighbour();

            scalarField magPhi = mag(phi);

            scalarField& wIn = w.internalField();

            // Larger small for complex arithmetic accuracy
            const scalar kSmall = 1000*SMALL;

            // Calculate internal weights using field magnitude
            forAll (owner, faceI)
            {
                label own = owner[faceI];
                label nei = neighbour[faceI];

                scalar mOwn = magPhi[own]/(1 - weights[faceI]);
                scalar mNei = magPhi[nei]/weights[faceI];

                scalar den = magPhi[nei] - magPhi[own];

                scalar mean = mOwn*mNei/
                    ((mOwn + mNei)*longDelta[faceI]*deltaCoeffs[faceI]);

                // Note: complex arithmetic requires extra accuracy
                // This is a division of two close subtractions
                // HJ, 28/Sep/2011
                if (mag(den) > kSmall)
                {
                    // mOwn > mNei
                    wIn[faceI] = (magPhi[nei] - mean)/den;
                }
                else
                {
                    wIn[faceI] = 0.5;
                }
            }

            forAll (phi.boundaryField(), pi)
            {
                fvsPatchScalarField& wp = w.boundaryField()[pi];

                const fvPatchField<Type>& pPhi = phi.boundaryField()[pi];

                if (pPhi.coupled())
                {
                    const scalarField& pWeights = weights.boundaryField()[pi];

                    scalarField magPhiOwn = mag(pPhi.patchInternalField());
                    scalarField magPhiNei = mag(pPhi.patchNeighbourField());

                    const scalarField& pDeltaCoeffs =
                        deltaCoeffs.boundaryField()[pi];

                    const scalarField& pLongDelta = mld.magDelta(pi);

                    // Calculate internal weights using field magnitude
                    forAll (pPhi, faceI)
                    {
                        scalar mOwn = magPhiOwn[faceI]/(1 - pWeights[faceI]);
                        scalar mNei = magPhiNei[faceI]/pWeights[faceI];

                        scalar den = magPhiNei[faceI] - magPhiOwn[faceI];

                        // Note: complex arithmetic requires extra accuracy
                        // This is a division of two close subtractions
                        // HJ, 28/Sep/2011
                        if (mag(den) > kSmall)
                        {
                            scalar mean = mOwn*mNei/
                                (
                                    (mOwn + mNei)*
                                    pLongDelta[faceI]*
                                    pDeltaCoeffs[faceI]
                                );

                            wp[faceI] = (magPhiNei[faceI] - mean)/den;
                        }
                        else
                        {
                            wp[faceI] = 0.5;
                        }
                    }
                }
                else
                {
                    // Boundary weights for uncoupled patches are 1
                    wp = 1;
                }
            }

            return tw;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
