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
    Foam::internalWriter

Description
    Write fields (internal).

SourceFiles
    internalWriter.C
    internalWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef internalWriter_H
#define internalWriter_H

#include "OFstream.H"
#include "volFields.H"
#include "pointFields.H"
#include "vtkMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class volPointInterpolation;

/*---------------------------------------------------------------------------*\
                           Class internalWriter Declaration
\*---------------------------------------------------------------------------*/

class internalWriter
{
    const vtkMesh& vMesh_;

    const bool binary_;

    const fileName fName_;

    std::ofstream os_;


public:

    // Constructors

        //- Construct from components
        internalWriter
        (
            const vtkMesh&,
            const bool binary,
            const fileName&
        );


    // Member Functions

        std::ofstream& os()
        {
            return os_;
        }

        //- Write cellIDs
        void writeCellIDs();

        //- Write generic GeometricFields
        template<class Type, template<class> class PatchField, class GeoMesh>
        void write
        (
            const PtrList<GeometricField<Type, PatchField, GeoMesh> >&
        );

        //- Interpolate and write volFields
        template<class Type>
        void write
        (
            const volPointInterpolation&,
            const PtrList<GeometricField<Type, fvPatchField, volMesh> >&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "internalWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
