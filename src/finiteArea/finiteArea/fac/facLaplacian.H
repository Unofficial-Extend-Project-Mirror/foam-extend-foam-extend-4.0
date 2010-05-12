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
    fac

Description
    Calculate the laplacian of the given field.

SourceFiles
    facLaplacian.C

\*---------------------------------------------------------------------------*/


#ifndef facLaplacian_H
#define facLaplacian_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fac functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const dimensionedScalar&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const dimensionedScalar&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const dimensionedScalar&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const dimensionedScalar&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const areaScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<areaScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const areaScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<areaScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const areaScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<areaScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const areaScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<areaScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const edgeScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<edgeScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
        const word&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const edgeScalarField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<edgeScalarField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );

    /*
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const areaTensorField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<areaTensorField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const areaTensorField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<areaTensorField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );


    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const edgeTensorField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<edgeTensorField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const edgeTensorField&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > laplacian
    (
        const tmp<edgeTensorField>&,
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );
    */
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "facLaplacian.C"
#else
#   ifdef xlC
#   pragma implementation("facLaplacian.C")
#   endif
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
