/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "facLnGrad.H"
#include "faMesh.H"
#include "lnGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::lnGradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().schemesDict().lnGradScheme(name)
    )().lnGrad(vf);
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf,
    const word& name
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > LnGrad
    (
        fac::lnGrad(tvf(), name)
    );
    tvf.clear();
    return LnGrad;
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
lnGrad
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fac::lnGrad(vf, "lnGrad(" + vf.name() + ')');
}


template<class Type>
tmp<GeometricField<Type, faePatchField, edgeMesh> >
lnGrad
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tvf
)
{
    tmp<GeometricField<Type, faePatchField, edgeMesh> > LnGrad
    (
        fac::lnGrad(tvf())
    );
    tvf.clear();
    return LnGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
