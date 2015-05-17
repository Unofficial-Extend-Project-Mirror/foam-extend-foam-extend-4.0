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

Description


\*---------------------------------------------------------------------------*/

#include "areaFields.H"
#include "edgeFields.H"
#include "faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type> >
Su
(
    const GeometricField<Type, faPatchField, areaMesh>& su,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*su.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam();

    fam.source() -= mesh.S()*su.internalField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type> >
Su
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& tsu,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > tfam = fam::Su(tsu(), vf);
    tsu.clear();
    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
Sp
(
    const areaScalarField& sp,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam();

    fam.diag() += mesh.S()*sp.internalField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type> >
Sp
(
    const tmp<areaScalarField>& tsp,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > tfam = fam::Sp(tsp(), vf);
    tsp.clear();
    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
Sp
(
    const dimensionedScalar& sp,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam();

    fam.diag() += mesh.S()*sp.value();

    return tfam;
}


template<class Type>
tmp<faMatrix<Type> >
SuSp
(
    const areaScalarField& sp,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    const faMesh& mesh = vf.mesh();

    tmp<faMatrix<Type> > tfam
    (
        new faMatrix<Type>
        (
            vf,
            dimArea*sp.dimensions()*vf.dimensions()
        )
    );
    faMatrix<Type>& fam = tfam();

    fam.diag() += mesh.S()*max(sp.internalField(), scalar(0));

    fam.source() -= mesh.S()*min(sp.internalField(), scalar(0))
        *vf.internalField();

    return tfam;
}

template<class Type>
tmp<faMatrix<Type> >
SuSp
(
    const tmp<areaScalarField>& tsp,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type> > tfam = fam::SuSp(tsp(), vf);
    tsp.clear();
    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
