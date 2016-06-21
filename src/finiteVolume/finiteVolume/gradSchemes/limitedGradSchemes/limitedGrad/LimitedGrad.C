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

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class GradientLimiter>
tmp<GeometricField<Type, fvPatchField, volMesh> >
LimitedGrad<Type, GradientLimiter>::limiter
(
    const GeoFieldType& vf,
    const GeoGradFieldType& gradVf
) const
{
    // Get reference to the mesh
    const fvMesh& mesh = vf.mesh();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    // Calculate min/max of field

    Field<Type> maxVf(vf.internalField());
    Field<Type> minVf(vf.internalField());

    const Field<Type>& vfIn = vf.internalField();

    // Internal faces
    forAll (owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        // Owner side
        maxVf[own] = Foam::max(maxVf[own], vfIn[nei]);
        minVf[own] = Foam::min(minVf[own], vfIn[nei]);

        // Neighbour side
        maxVf[nei] = Foam::max(maxVf[nei], vfIn[own]);
        minVf[nei] = Foam::min(minVf[nei], vfIn[own]);
    }

    const GeoBoundaryFieldType& bf = vf.boundaryField();

    // Boundary faces
    forAll (bf, patchI)
    {
        const fvPatchField<Type>& psf = bf[patchI];

        const unallocLabelList& pOwner = mesh.boundary()[patchI].faceCells();

        if (psf.coupled())
        {
            // For a coupled boundary, use neighbour field
            Field<Type> psfNei = psf.patchNeighbourField();

            forAll (pOwner, pFaceI)
            {
                label own = pOwner[pFaceI];

                maxVf[own] = Foam::max(maxVf[own], psfNei[pFaceI]);
                minVf[own] = Foam::min(minVf[own], psfNei[pFaceI]);
            }
        }
        else
        {
            // For regular boundary, use boundary value
            forAll (pOwner, pFaceI)
            {
                label own = pOwner[pFaceI];
                const Type& vfNei = psf[pFaceI];

                maxVf[own] = Foam::max(maxVf[own], vfNei);
                minVf[own] = Foam::min(minVf[own], vfNei);
            }
        }
    }

    // Subtract the cell value to get differences
    // Stabilise differences for round-off error, since maxVf must be
    // positive or zero and minVf must be negative or zero
    // HJ, 1/Feb/2016
    maxVf = Foam::max(maxVf - vf, pTraits<Type>::zero);
    minVf = Foam::min(minVf - vf, pTraits<Type>::zero);

    // Create a limiter
    tmp<GeoFieldType> tlimitField
    (
        new GeoFieldType
        (
            IOobject
            (
                "limitField(" + vf.name() + ")",
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensioned<Type>("zero", dimless, pTraits<Type>::one),
            zeroGradientFvPatchField<Type>::typeName
        )
    );
    GeoFieldType& limitField = tlimitField();

    const volVectorField& C = mesh.C();
    const vectorField& CIn = C.internalField();

    const surfaceVectorField& Cf = mesh.Cf();
    const vectorField& CfIn = Cf.internalField();

    const scalarField& vols = mesh.V().field();

    Field<Type>& lfIn = limitField.internalField();

    const GradFieldType& g = gradVf.internalField();

    // Apply limiter function
    GradientLimiter limitFunction;

    forAll (owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        // Owner side
        limitFunction.limiter
        (
            lfIn[own],
            vols[own],
            maxVf[own],
            minVf[own],
            (CfIn[faceI] - CIn[own]) & g[own]
        );

        // Neighbour side
        limitFunction.limiter
        (
            lfIn[nei],
            vols[nei],
            maxVf[nei],
            minVf[nei],
            (CfIn[faceI] - CIn[nei]) & g[nei]
        );
    }

    forAll (bf, patchI)
    {
        const unallocLabelList& pOwner = mesh.boundary()[patchI].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchI];

        forAll (pOwner, pFaceI)
        {
            label own = pOwner[pFaceI];

            limitFunction.limiter
            (
                lfIn[own],
                vols[own],
                maxVf[own],
                minVf[own],
                (pCf[pFaceI] - C[own]) & g[own]
            );
        }
    }

    // Update coupled boundaries for patchNeighbourField
    limitField.correctBoundaryConditions();

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vf.name()
            << " max = " << gMax(lfIn)
            << " min = " << gMin(lfIn)
            << " average: " << gAverage(lfIn) << endl;
    }

    return tlimitField;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GradientLimiter>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
LimitedGrad<Type, GradientLimiter>::gradientField
(
    const GeoFieldType& vf,
    const word& name
) const
{
    // Get base gradient
    tmp<GeoGradFieldType> tGrad = basicGradScheme_().calcGrad(vf, name);
    GeoGradFieldType& gradVf = tGrad();

    // Apply the limiter
    GeoFieldType limitField(this->limiter(vf, gradVf));

    gradVf.internalField() =
        scaleRow(gradVf.internalField(), limitField.internalField());

    gradVf.correctBoundaryConditions();
    gradScheme<Type>::correctBoundaryConditions(vf, gradVf);

    return tGrad;
}


template<class Type, class GradientLimiter>
tmp
<
    BlockLduSystem<vector, typename outerProduct<vector, Type>::type>
>
LimitedGrad<Type, GradientLimiter>::gradientMatrix
(
   const GeoFieldType& vf
) const
{
    // Calculate base gradient matrix
    tmp<GradMatrixType> tbs = basicGradScheme_().fvmGrad(vf);
    GradMatrixType& bs = tbs();

    // Calculate limiter.  Using explicit gradient
    // Using cached gradient?  Check.  HJ, 4/Jun/2016
    GeoFieldType limitField
    (
        this->limiter
        (
            vf,
            basicGradScheme_().grad(vf)()
        )
    );

    const Field<Type>& lfIn = limitField.internalField();

    typedef typename CoeffField<vector>::linearTypeField
        linearCoeffType;

    typedef typename outerProduct<vector, Type>::type sourceType;

    Field<sourceType>& source = bs.source();

    // Grab ldu parts of block matrix as linear always
    linearCoeffType& d = bs.diag().asLinear();

    linearCoeffType& u = bs.upper().asLinear();

    linearCoeffType& l = bs.lower().asLinear();

    // Limit upper and lower coeffs

    const fvMesh& mesh = vf.mesh();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll (u, faceI)
    {
        u[faceI] = scaleRow(u[faceI], lfIn[owner[faceI]]);
        l[faceI] = scaleRow(l[faceI], lfIn[neighbour[faceI]]);
    }

    // Limit diag and source coeffs
    forAll (d, cellI)
    {
        d[cellI] = scaleRow(d[cellI], lfIn[cellI]);

        source[cellI] = scaleRow(source[cellI], lfIn[cellI]);
    }

    // Limit coupling coeffs
    forAll (vf.boundaryField(), patchI)
    {
        const fvPatchField<Type>& pf = vf.boundaryField()[patchI];
        const fvPatch& patch = pf.patch();

        const labelList& fc = patch.faceCells();

        if (patch.coupled())
        {
            linearCoeffType& pcoupleUpper =
                bs.coupleUpper()[patchI].asLinear();

            linearCoeffType& pcoupleLower =
                bs.coupleLower()[patchI].asLinear();

            const Field<Type> lfNei =
                limitField.boundaryField()[patchI].patchNeighbourField();

            forAll (pf, faceI)
            {
                pcoupleUpper[faceI] =
                    scaleRow(pcoupleUpper[faceI], lfIn[fc[faceI]]);

                pcoupleLower[faceI] =
                    scaleRow(pcoupleLower[faceI], lfNei[faceI]);
            }
        }
    }

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
