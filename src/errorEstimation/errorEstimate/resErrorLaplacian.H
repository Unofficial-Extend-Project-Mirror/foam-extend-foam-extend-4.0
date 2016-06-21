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
    Foam::resError

Description
    Residual error estimate for the fv laplacian operators

SourceFiles
    resErrorLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef resErrorLaplacian_H
#define resErrorLaplacian_H

#include "errorEstimate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace resError
{
    // Laplacian terms

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const dimensionedScalar&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const volScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const tmp<volScalarField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const surfaceScalarField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const tmp<surfaceScalarField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const volTensorField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const tmp<volTensorField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const surfaceTensorField&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<errorEstimate<Type> > laplacian
        (
            const tmp<surfaceTensorField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "resErrorLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
