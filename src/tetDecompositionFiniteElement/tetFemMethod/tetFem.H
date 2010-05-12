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

Class
    tetFem

Description
    Class of static functions to calculate implicit finite element derivatives
    returning a matrix.

SourceFiles
    tetFem.C

\*---------------------------------------------------------------------------*/

#ifndef tetFem_H
#define tetFem_H

#include "tetPointFieldsFwd.H"
#include "tetFemMatrices.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Namespace tetFem functions Declaration
\*---------------------------------------------------------------------------*/

namespace tetFem
{
    // Laplacian

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacian
        (
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacian
        (
            const elementScalarField&,
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > smoother
        (
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacian
        (
            const dimensionedScalar&,
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacianTranspose
        (
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacianTranspose
        (
            const elementScalarField&,
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacianTranspose
        (
            const dimensionedScalar&,
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacianTrace
        (
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );
        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacianTrace
        (
            const elementScalarField&,
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );

        template<class Type>
        static tmp<tetFemMatrix<Type> > laplacianTrace
        (
            const dimensionedScalar&,
            GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetFem.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
