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
    Edge integrate edgeField creating a areaField.
    Edge sum a edgeField creating a areaField.

SourceFiles
    facEdgeIntegrate.C

\*---------------------------------------------------------------------------*/


#ifndef facEdgeIntegrate_H
#define facEdgeIntegrate_H

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
    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> >
    edgeIntegrate
    (
        const GeometricField<Type, faePatchField, edgeMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> >
    edgeIntegrate
    (
        const tmp<GeometricField<Type, faePatchField, edgeMesh> >&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > edgeSum
    (
        const GeometricField<Type, faePatchField, edgeMesh>&
    );

    template<class Type>
    tmp<GeometricField<Type, faPatchField, areaMesh> > edgeSum
    (
        const tmp<GeometricField<Type, faePatchField, edgeMesh> >&
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "facEdgeIntegrate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
