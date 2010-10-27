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

#include "gaussFaGrad.H"
#include "facGrad.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, faPatchField, areaMesh
    >
>
gaussGrad<Type>::grad
(
    const GeometricField<Type, faPatchField, areaMesh>& vsf
) const
{
    typedef typename outerProduct<vector, Type>::type GradType;

    tmp<GeometricField<GradType, faPatchField, areaMesh> > tgGrad
    (
        fac::edgeIntegrate
        (
            vsf.mesh().Le()
           *tinterpScheme_().interpolate(vsf)
        )
    );

    GeometricField<GradType, faPatchField, areaMesh>& gGrad = tgGrad();

    // Remove component of gradient normal to surface (area)
    const areaVectorField& n = vsf.mesh().faceAreaNormals();
//     const areaScalarField& K = vsf.mesh().faceCurvatures();

    gGrad -= n*(n & gGrad);
//     gGrad -= K*(n*vsf);
    gGrad.correctBoundaryConditions();

    gGrad.rename("grad(" + vsf.name() + ')');
    correctBoundaryConditions(vsf, gGrad);

    return tgGrad;
}


template<class Type>
void gaussGrad<Type>::correctBoundaryConditions
(
    const GeometricField<Type, faPatchField, areaMesh>& vsf,
    GeometricField
    <
        typename outerProduct<vector, Type>::type, faPatchField, areaMesh
    >& gGrad
)
{
    forAll(vsf.boundaryField(), patchI)
    {
        if (!vsf.boundaryField()[patchI].coupled())
        {
            vectorField m =
                vsf.mesh().Le().boundaryField()[patchI]
                /vsf.mesh().magLe().boundaryField()[patchI];

            // Zeljko Tukovic
//             gGrad.boundaryField()[patchI] = 
//                 m*vsf.boundaryField()[patchI].snGrad();

            //HJ Not sure: should this be a complete correction or just the
            //   tangential part?  HJ, 24/Jul/2009
            gGrad.boundaryField()[patchI] += m*
            (
                vsf.boundaryField()[patchI].snGrad()
              - (m & gGrad.boundaryField()[patchI])
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
