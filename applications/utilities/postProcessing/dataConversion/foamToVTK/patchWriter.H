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
    Foam::patchWriter

Description
    Write patch fields

SourceFiles
    patchWriter.C
    patchWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef patchWriter_H
#define patchWriter_H

#include "pointMesh.H"
#include "OFstream.H"
#include "volFields.H"
#include "pointFields.H"
#include "vtkMesh.H"
#include "indirectPrimitivePatch.H"
#include "PrimitivePatchInterpolationTemplate.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class volPointInterpolation;

/*---------------------------------------------------------------------------*\
                           Class patchWriter Declaration
\*---------------------------------------------------------------------------*/

class patchWriter
{
    const vtkMesh& vMesh_;

    const bool binary_;

    const bool nearCellValue_;

    const fileName fName_;

    const labelList patchIDs_;

    std::ofstream os_;

    label nPoints_;

    label nFaces_;


public:

    // Constructors

        //- Construct from components
        patchWriter
        (
            const vtkMesh&,
            const bool binary,
            const bool nearCellValue,
            const fileName&,
            const labelList& patchIDs
        );


    // Member Functions

        std::ofstream& os()
        {
            return os_;
        }

        label nPoints() const
        {
            return nPoints_;
        }

        label nFaces() const
        {
            return nFaces_;
        }

        //- Write cellIDs
        void writePatchIDs();

        //- Write volFields
        template<class Type>
        void write
        (
            const PtrList<GeometricField<Type, fvPatchField, volMesh> >&
        );

        //- Write pointFields
        template<class Type>
        void write
        (
            const PtrList<GeometricField<Type, pointPatchField, pointMesh> >&
        );

        //- Interpolate and write volFields
        template<class Type>
        void write
        (
            const PrimitivePatchInterpolation<primitivePatch>&,
            const PtrList<GeometricField<Type, fvPatchField, volMesh> >&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "patchWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
