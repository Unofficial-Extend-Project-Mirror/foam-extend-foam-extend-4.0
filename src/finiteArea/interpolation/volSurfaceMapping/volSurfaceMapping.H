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
    volSurfaceMapping

Description
    Volume to surface and surface to volume mapping

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    volSurfaceMapping.C

\*---------------------------------------------------------------------------*/

#ifndef volSurfaceMapping_H
#define volSurfaceMapping_H

#include "faMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvPatchField;

/*---------------------------------------------------------------------------*\
                        Class volSurfaceMapping Declaration
\*---------------------------------------------------------------------------*/

class volSurfaceMapping
{
    // Private data

        //- Reference to mesh
        const faMesh& mesh_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        volSurfaceMapping(const volSurfaceMapping&);

        //- Disallow default bitwise assignment
        void operator=(const volSurfaceMapping&);


public:

    // Constructors

        //- Construct null
        volSurfaceMapping();

        //- Construct from mesh
        volSurfaceMapping(const faMesh& mesh)
        :
            mesh_(mesh)
        {}


    // Destructor - default


    // Member Functions

        //- Map droplet cloud sources to surface
        template<template<class> class PatchField, class Type>
        tmp<Field<Type> > mapToSurface
        (
            const FieldField<PatchField, Type>& df
        ) const;

        //- Map surface field to volume boundary field
        template<class Type>
        void mapToVolume
        (
            const GeometricField<Type, faPatchField, areaMesh>& af,
            FieldField<fvPatchField, Type>& bf
        ) const;

        template<class Type>
        void mapToVolume
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh> >& taf,
            FieldField<fvPatchField, Type>& bf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "volSurfaceMapping.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
