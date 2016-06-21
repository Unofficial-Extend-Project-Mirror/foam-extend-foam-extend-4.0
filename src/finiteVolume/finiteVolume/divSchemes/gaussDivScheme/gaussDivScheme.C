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

#include "gaussDivScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
>
gaussDivScheme<Type>::fvcDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp
    <
        GeometricField
        <typename innerProduct<vector, Type>::type, fvPatchField, volMesh>
    > tDiv
    (
        fvc::surfaceIntegrate
        (
            this->mesh_.Sf() & this->tinterpScheme_().interpolate(vf)
        )
    );

    tDiv().rename("div(" + vf.name() + ')');

    return tDiv;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> gaussDivScheme<Type>::fvmUDiv
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> gaussDivScheme<Type>::fvmUDiv\n"
        "(\n"
        "    GeometricField<Type, fvPatchField, volMesh>&"
        ")\n"
    )   << "Implicit div operator defined only for vector."
        << abort(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType> > tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


template<class Type>
tmp
<
    BlockLduSystem<vector, typename innerProduct<vector, Type>::type>
> gaussDivScheme<Type>::fvmUDiv
(
    const surfaceScalarField& flux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    FatalErrorIn
    (
        "tmp<BlockLduSystem> gaussDivScheme<Type>::fvmUDiv\n"
        "(\n"
        "    const surfaceScalarField& flux"
        "    const GeometricField<Type, fvPatchField, volMesh>&"
        ")\n"
    )   << "Implicit div operator defined only for vector."
        << abort(FatalError);

    typedef typename innerProduct<vector, Type>::type DivType;

    tmp<BlockLduSystem<vector, DivType> > tbs
    (
        new BlockLduSystem<vector, DivType>(vf.mesh())
    );

    return tbs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
