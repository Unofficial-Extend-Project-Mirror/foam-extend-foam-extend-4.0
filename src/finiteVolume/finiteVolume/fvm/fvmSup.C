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

#include "fvmSup.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::Su
(
    const DimensionedField<Type, volMesh>& su,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*su.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    fvm.source() -= mesh.V()*su.field();

    return tfvm;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::Su
(
    const tmp<DimensionedField<Type, volMesh> >& tsu,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::Su
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = fvm::Su(tsu(), vf);
    tsu.clear();
    return tfvm;
}

template<class Type>
Foam::zeroField
Foam::fvm::Su
(
    const zeroField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::Sp
(
    const DimensionedField<scalar, volMesh>& sp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    fvm.diag() += mesh.V()*sp.field();

    return tfvm;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::Sp
(
    const tmp<DimensionedField<scalar, volMesh> >& tsp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::Sp
(
    const tmp<volScalarField>& tsp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = fvm::Sp(tsp(), vf);
    tsp.clear();
    return tfvm;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::Sp
(
    const dimensionedScalar& sp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*sp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    fvm.diag() += mesh.V()*sp.value();

    return tfvm;
}

template<class Type>
Foam::zeroField
Foam::fvm::Sp
(
    const zeroField&,
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return zeroField();
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::SuSp
(
    const DimensionedField<scalar, volMesh>& susp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const fvMesh& mesh = vf.mesh();

    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            dimVol*susp.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm();

    fvm.diag() += mesh.V()*max(susp.field(), scalar(0));

    fvm.source() -= mesh.V()*min(susp.field(), scalar(0))
        *vf.internalField();

    return tfvm;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::SuSp
(
    const tmp<DimensionedField<scalar, volMesh> >& tsusp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}

template<class Type>
Foam::tmp<Foam::fvMatrix<Type> >
Foam::fvm::SuSp
(
    const tmp<volScalarField>& tsusp,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type> > tfvm = fvm::SuSp(tsusp(), vf);
    tsusp.clear();
    return tfvm;
}

template<class Type>
Foam::zeroField
Foam::fvm::SuSp
(
    const zeroField&,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return zeroField();
}


// ************************************************************************* //
