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

Description
    
\*---------------------------------------------------------------------------*/

#include "gaussFaLaplacianScheme.H"
#include "facDiv.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type> >
gaussLaplacianScheme<Type>::famLaplacian
(
    const edgeScalarField& gamma,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<edgeScalarField> tdeltaCoeffs = this->tlnGradScheme_().deltaCoeffs(vf);
    const edgeScalarField& deltaCoeffs = tdeltaCoeffs();

    edgeScalarField gammaMagSf = gamma*this->mesh().magLe();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            deltaCoeffs.dimensions()*gammaMagSf.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam();

    fam.upper() = deltaCoeffs.internalField()*gammaMagSf.internalField();
    fam.negSumDiag();

    forAll(fam.psi().boundaryField(), patchI)
    {
        const faPatchField<Type>& psf = fam.psi().boundaryField()[patchI];
        const faePatchScalarField& patchGamma =
            gammaMagSf.boundaryField()[patchI];

        fam.internalCoeffs()[patchI] = patchGamma*psf.gradientInternalCoeffs();
        fam.boundaryCoeffs()[patchI] = -patchGamma*psf.gradientBoundaryCoeffs();
    }

    if (this->tlnGradScheme_().corrected())
    {
        if (this->mesh().fluxRequired(vf.name()))
        {
            fam.faceFluxCorrectionPtr() = new
            GeometricField<Type, faePatchField, edgeMesh>
            (
                gammaMagSf*this->tlnGradScheme_().correction(vf)
            );

            fam.source() -=
                this->mesh().S()*
                fac::div
                (
                    *fam.faceFluxCorrectionPtr()
                )().internalField();
        }
        else
        {
            fam.source() -=
                this->mesh().S()*
                fac::div
                (
                    gammaMagSf*this->tlnGradScheme_().correction(vf)
                )().internalField();
        }
    }

    return tfam;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
gaussLaplacianScheme<Type>::facLaplacian
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > tLaplacian
    (
        fac::div(this->tlnGradScheme_().lnGrad(vf)*vf.mesh().magLe())
    );

    tLaplacian().rename("laplacian(" + vf.name() + ')');

    return tLaplacian;
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh> >
gaussLaplacianScheme<Type>::facLaplacian
(
    const edgeScalarField& gamma,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<GeometricField<Type, faPatchField, areaMesh> > tLaplacian
    (
        fac::div(gamma*this->tlnGradScheme_().lnGrad(vf)*vf.mesh().magLe())
    );

    tLaplacian().rename("laplacian(" + gamma.name() + ',' + vf.name() + ')');

    return tLaplacian;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
