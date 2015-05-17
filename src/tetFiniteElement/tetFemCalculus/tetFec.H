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

Class
    tetFec

Description
    Class of static functions to calculate explicit finite element derivatives.

SourceFiles
    tetFec.C

\*---------------------------------------------------------------------------*/

#ifndef tetFec_H
#define tetFec_H

#include "tetPointFieldsFwd.H"
#include "elementFieldsFwd.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Namespace tetFec functions Declaration
\*---------------------------------------------------------------------------*/

namespace tetFec
{
    // Gradient

        template<class Type>
        static tmp
        <
            GeometricField
            <
                typename outerProduct<vector, Type>::type,
                tetPolyPatchField,
                tetPointMesh
            >
        > grad
        (
            const GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );


        template<class Type>
        static tmp
        <
            GeometricField
            <
                typename outerProduct<vector, Type>::type,
                elementPatchField,
                elementMesh
            >
        > elementGrad
        (
            const GeometricField<Type, tetPolyPatchField, tetPointMesh>&
        );


        template<class Type>
        tmp<GeometricField<Type, elementPatchField, elementMesh> > ddt
        (
            const GeometricField<Type, elementPatchField, elementMesh>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetFec.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
