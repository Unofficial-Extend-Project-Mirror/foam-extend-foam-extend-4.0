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

InNamespace
    Foam::fac

Description
    Edge-weighted average a edgeField creating a areaField

SourceFiles
    facAverage.C

\*---------------------------------------------------------------------------*/


#ifndef facAverage_H
#define facAverage_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fac functions Declaration
\*---------------------------------------------------------------------------*/

namespace fac
{
    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > average
    (
        const GeometricField<Type,  faePatchField, edgeMesh>&
    );

    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > average
    (
        const tmp<GeometricField<Type,  faePatchField, edgeMesh> >&
    );


    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > average
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    //- Area-weighted average a edgeField creating a areaField
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > average
    (
        const tmp<GeometricField<Type, faPatchField, areaMesh> >&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "facAverage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
