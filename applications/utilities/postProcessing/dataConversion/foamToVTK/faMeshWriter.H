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
#include "PrimitivePatchInterpolationTemplate.H"

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
