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

InNamespace
    Foam::fvm

Description
    Calculate the matrix for implicit and explicit sources.

SourceFiles
    fvmSup.C

\*---------------------------------------------------------------------------*/

#ifndef fvmSup_H
#define fvmSup_H

#include "volFieldsFwd.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvm functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvm
{
   // Explicit source

        template<class Type>
        tmp<fvMatrix<Type> > Su
        (
            const DimensionedField<Type, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type> > Su
        (
            const tmp<DimensionedField<Type, volMesh> >&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type> > Su
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        zeroField Su
        (
            const zeroField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // Implicit source

        template<class Type>
        tmp<fvMatrix<Type> > Sp
        (
            const DimensionedField<scalar, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type> > Sp
        (
            const tmp<DimensionedField<scalar, volMesh> >&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type> > Sp
        (
            const tmp<volScalarField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        template<class Type>
        tmp<fvMatrix<Type> > Sp
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


        template<class Type>
        zeroField Sp
        (
            const zeroField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        template<class Type>
        tmp<fvMatrix<Type> > SuSp
        (
            const DimensionedField<scalar, volMesh>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type> > SuSp
        (
            const tmp<DimensionedField<scalar, volMesh> >&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type> > SuSp
        (
            const tmp<volScalarField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        zeroField SuSp
        (
            const zeroField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvmSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
