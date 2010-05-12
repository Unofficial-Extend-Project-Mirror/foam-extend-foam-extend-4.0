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
    Foam::surfaceWriter

Description

SourceFiles
    surfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWriter_H
#define surfaceWriter_H

#include "Field.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class surfaceWriter
{

public:

    //- Runtime type information
    TypeName("surfaceWriter");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            word,
            (),
            ()
        );


    // Selectors

        //- Return a reference to the selected surfaceWriter
        static autoPtr<surfaceWriter> New(const word& writeFormat);


    // Constructors

        //- Construct null
        surfaceWriter();


    // Destructor

        virtual ~surfaceWriter();


    // Member Functions

        //- Writes single surface to file. Either one value per vertex or
        //  one value per face (detected by values.size()==faces.size())
        virtual void write
        (
            const fileName& samplePath,     // <root>/<case>/sampleSurfaces
            const fileName& timeDir,        // time directory
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const fileName& fieldName,      // name of field
            const Field<Type>& values,
            const bool verbose = false
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "surfaceWriter.C"
#endif

// Only used internally
#define makeTypeSurfaceWritersTypeName(type)                                 \
                                                                             \
defineNamedTemplateTypeNameAndDebug(type, 0);

// Used externally sometimes
#define makeSurfaceWritersTypeName(typeWriter)                               \
                                                                             \
makeTypeSurfaceWritersTypeName(typeWriter##ScalarWriter);                    \
makeTypeSurfaceWritersTypeName(typeWriter##VectorWriter);                    \
makeTypeSurfaceWritersTypeName(typeWriter##SphericalTensorWriter);           \
makeTypeSurfaceWritersTypeName(typeWriter##SymmTensorWriter);                \
makeTypeSurfaceWritersTypeName(typeWriter##TensorWriter);

// Define type info for single template instantiation (e.g. vector)
#define makeSurfaceWriterTypes(WriterType, type)                             \
                                                                             \
defineNamedTemplateTypeNameAndDebug(type, 0);                                \
                                                                             \
addToRunTimeSelectionTable                                                   \
(                                                                            \
    WriterType, type, word                                                   \
);


// Define type info info for scalar, vector etc. instantiations
#define makeSurfaceWriters(typeWriter)                                       \
                                                                             \
makeSurfaceWriterTypes(scalarWriter, typeWriter##ScalarWriter);              \
makeSurfaceWriterTypes(vectorWriter, typeWriter##VectorWriter);              \
makeSurfaceWriterTypes(sphericalTensorWriter, typeWriter##SphericalTensorWriter);\
makeSurfaceWriterTypes(symmTensorWriter, typeWriter##SymmTensorWriter);      \
makeSurfaceWriterTypes(tensorWriter, typeWriter##TensorWriter);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
