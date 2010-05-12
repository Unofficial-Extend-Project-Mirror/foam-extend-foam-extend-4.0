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
    Fourth-order snGrad scheme with non-orthogonal correction.

\*---------------------------------------------------------------------------*/

#include "fourthLnGrad.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "correctedLnGrad.H"
#include "linearEdgeInterpolation.H"
#include "gaussFaGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
fourthLnGrad<Type>::~fourthLnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
fourthLnGrad<Type>::correction
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
) const
{
    const faMesh& mesh = this->mesh();

    tmp<GeometricField<Type, faePatchField, edgeMesh> > tcorr
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                "lnGradCorr("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()*this->mesh().deltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, faePatchField, edgeMesh>& corr = tcorr();

    edgeVectorField m = mesh.Le()/mesh.magLe();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; cmpt++)
    {
        corr.replace
        (
            cmpt,
          - (1.0/15.0)*m
          & linearEdgeInterpolation
            <
                typename
                outerProduct<vector, typename pTraits<Type>::cmptType>::type
            >(mesh).interpolate
            (
                gaussGrad<typename pTraits<Type>::cmptType>(mesh)
               .grad(vf.component(cmpt))
            )
        );
    }

    corr += (1.0/15.0)*correctedLnGrad<Type>(mesh).lnGrad(vf);


//     tmp<GeometricField<Type, faePatchField, edgeMesh> > tcorr
//     (
//         (1.0/15.0)
//        *(
//             correctedLnGrad<Type>(mesh).lnGrad(vf)
//           - (
//               linearEdgeInterpolate(gaussGrad<Type>(mesh).grad(vf)) & mesh.Le()
//             )/mesh.magLe()
//         )
//     );

    if (correctedLnGrad<Type>(mesh).corrected())
    {
        tcorr() += correctedLnGrad<Type>(mesh).correction(vf);
    }

    return tcorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
