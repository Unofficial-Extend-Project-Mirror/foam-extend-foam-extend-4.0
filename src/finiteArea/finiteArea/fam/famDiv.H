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
    fam

Description
    Calculate the matrix for the divergence of the given field and flux.

SourceFiles
    famDiv.C

\*---------------------------------------------------------------------------*/

#ifndef famDiv_H
#define famDiv_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "edgeInterpolationScheme.H"
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
    tmp<faMatrix<Type> > div
    (
        const edgeScalarField&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<faMatrix<Type> > div
    (
        const tmp<edgeScalarField>&,
        GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<faMatrix<Type> > div
    (
        const edgeScalarField&,
        GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type> > div
    (
        const tmp<edgeScalarField>&,
        GeometricField<Type, faPatchField, areaMesh>&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "famDiv.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
