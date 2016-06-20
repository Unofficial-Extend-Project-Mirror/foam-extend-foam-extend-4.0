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

Description
    Residual error estimate for the fv laplacian operators.

\*---------------------------------------------------------------------------*/

#include "resErrorLaplacian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace resError
{

template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            "gamma",
            vf.time().constant(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        dimensionedScalar("1", dimless, 1.0)
    );

    return resError::laplacian(Gamma, vf);
}


template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const dimensionedScalar& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    surfaceScalarField Gamma
    (
        IOobject
        (
            gamma.name(),
            vf.time().timeName(),
            vf.db(),
            IOobject::NO_READ
        ),
        vf.mesh(),
        gamma
    );

    return resError::laplacian(Gamma, vf);
}


template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const volScalarField& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return resError::laplacian(fvc::interpolate(gamma), vf);
}


template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const tmp<volScalarField>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<errorEstimate<Type> > Laplacian(resError::laplacian(tgamma(), vf));
    tgamma.clear();
    return Laplacian;
}


template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const surfaceScalarField& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    const scalarField& vols = mesh.V();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField magSf = mesh.magSf();
    const fvPatchList& patches = mesh.boundary();
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    const surfaceScalarField& delta =
        mesh.surfaceInterpolation::deltaCoeffs();

    Field<Type> res(vols.size(), pTraits<Type>::zero);
    scalarField aNorm(vols.size(), 0.0);

    // Calculate gradient of the solution.
    // Change of return type due to gradient cacheing.  HJ, 22/Apr/2016
    const tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type,
            fvPatchField,
            volMesh
        >
    > tgradVf = fvc::grad(vf);

    const GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >& gradVf = tgradVf();

    // Internal faces
    forAll (owner, faceI)
    {
        // Owner

        // Subtract diffusion
        res[owner[faceI]] -=
            gamma[faceI]*(Sf[faceI] & gradVf[owner[faceI]]);

        aNorm[owner[faceI]] += delta[faceI]*gamma[faceI]*magSf[faceI];

        // Neighbour

        // Subtract diffusion
        res[neighbour[faceI]] +=
            gamma[faceI]*(Sf[faceI] & gradVf[neighbour[faceI]]);

        aNorm[neighbour[faceI]] += delta[faceI]*gamma[faceI]*magSf[faceI];

    }

    forAll (patches, patchI)
    {
        const vectorField& patchSf = Sf.boundaryField()[patchI];
        const scalarField& patchMagSf = magSf.boundaryField()[patchI];
        const scalarField& patchGamma = gamma.boundaryField()[patchI];
        const scalarField& patchDelta = delta.boundaryField()[patchI];

        const labelList& fCells = patches[patchI].faceCells();

        forAll (fCells, faceI)
        {
            // Subtract diffusion
            res[fCells[faceI]] -=
                patchGamma[faceI]*
                (
                    patchSf[faceI] & gradVf[fCells[faceI]]
                );

            aNorm[fCells[faceI]] +=
                patchDelta[faceI]*patchGamma[faceI]*patchMagSf[faceI];
        }
    }

    res /= vols;
    aNorm /= vols;

    return tmp<errorEstimate<Type> >
    (
        new errorEstimate<Type>
        (
            vf,
            delta.dimensions()*gamma.dimensions()*magSf.dimensions()
            *vf.dimensions(),
            res,
            aNorm
        )
    );
}

template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const tmp<surfaceScalarField>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<errorEstimate<Type> > tresError(resError::laplacian(tgamma(), vf));
    tgamma.clear();
    return tresError;
}


template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const volTensorField& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    return resError::laplacian
    (
        (mesh.Sf() & fvc::interpolate(gamma) & mesh.Sf())
        /sqr(mesh.magSf()),
        vf
    );
}

template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const tmp<volTensorField>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<errorEstimate<Type> > Laplacian = resError::laplacian(tgamma(), vf);
    tgamma.clear();
    return Laplacian;
}


template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const surfaceTensorField& gamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    return resError::laplacian
    (
        (mesh.Sf() & gamma & mesh.Sf())/sqr(mesh.magSf()),
        vf
    );
}

template<class Type>
tmp<errorEstimate<Type> >
laplacian
(
    const tmp<surfaceTensorField>& tgamma,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<errorEstimate<Type> > Laplacian = resError::laplacian(tgamma(), vf);
    tgamma.clear();
    return Laplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace resError

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

