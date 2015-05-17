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

Namespace
    fam

Description
    Calculate the matrix for implicit and explicit sources.

SourceFiles
    famSup.C

\*---------------------------------------------------------------------------*/

#ifndef famSup_H
#define famSup_H

#include "areaFieldsFwd.H"
#include "faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fam functions Declaration
\*---------------------------------------------------------------------------*/

namespace fam
{
   // Explicit source

        template<class Type>
        tmp<faMatrix<Type> > Su
        (
            const GeometricField<Type, faPatchField, areaMesh>&,
            GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type> > Su
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh> >&,
            GeometricField<Type, faPatchField, areaMesh>&
        );


    // Implicit source

        template<class Type>
        tmp<faMatrix<Type> > Sp
        (
            const areaScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type> > Sp
        (
            const tmp<areaScalarField>&,
            GeometricField<Type, faPatchField, areaMesh>&
        );


        template<class Type>
        tmp<faMatrix<Type> > Sp
        (
            const dimensionedScalar&,
            GeometricField<Type, faPatchField, areaMesh>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        template<class Type>
        tmp<faMatrix<Type> > SuSp
        (
            const areaScalarField&,
            GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type> > SuSp
        (
            const tmp<areaScalarField>&,
            GeometricField<Type, faPatchField, areaMesh>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "famSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
