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

Namespace
    famLaplacian

Description
    Calculate the matrix for the laplacian of the field.

SourceFiles
    famLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef famLaplacian_H
#define famLaplacian_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fam functions Declaration
\*---------------------------------------------------------------------------*/

namespace fam
{
    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const dimensionedScalar&,
        GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const dimensionedScalar&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const areaScalarField&,
        GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const areaScalarField&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const areaScalarField&,
        GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const tmp<areaScalarField>&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );
    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const tmp<areaScalarField>&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const edgeScalarField&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const tmp<edgeScalarField>&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const edgeScalarField&,
        GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const tmp<edgeScalarField>&,
        GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const areaTensorField&,
        GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const tmp<areaTensorField>&,
        GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const edgeTensorField&,
        GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type> > laplacian
    (
        const tmp<edgeTensorField>&,
        GeometricField<Type, faPatchField, areaMesh>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "famLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
