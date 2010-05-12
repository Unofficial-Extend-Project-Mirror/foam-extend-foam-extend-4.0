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
    Foam::faMeshWriter

Description
    Write Finite Area fields

SourceFiles
    faMeshWriter.C
    faMeshWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef faMeshWriter_H
#define faMeshWriter_H

#include "faMesh.H"
#include "OFstream.H"
#include "areaFields.H"
#include "indirectPrimitivePatch.H"
#include "PrimitivePatchInterpolation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class volPointInterpolation;

/*---------------------------------------------------------------------------*\
                           Class faMeshWriter Declaration
\*---------------------------------------------------------------------------*/

class faMeshWriter
{
    //- Reference to FA mesh
    const faMesh& aMesh_;

    //- Write in binary format
    const bool binary_;

    //- Output stream
    std::ofstream os_;


public:

    // Constructors

        //- Construct from components
        faMeshWriter
        (
            const faMesh&,
            const bool binary,
            const fileName&
        );


    // Member Functions

        //- Return stream
        std::ofstream& os()
        {
            return os_;
        }

        //- Write areaField
        template<class Type>
        void write
        (
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        //- Write list of area fields
        template<class Type>
        void write
        (
            const PtrList<GeometricField<Type, faPatchField, areaMesh> >&
        );

        //- Interpolate and write areaField
        template<class Type>
        void write
        (
            const PrimitivePatchInterpolation<indirectPrimitivePatch>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        //- Interpolate and write areaFields
        template<class Type>
        void write
        (
            const PrimitivePatchInterpolation<indirectPrimitivePatch>&,
            const PtrList<GeometricField<Type, faPatchField, areaMesh> >&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "faMeshWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
