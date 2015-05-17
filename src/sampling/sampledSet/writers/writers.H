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

InClass
    Foam::writer

SourceFiles
    writers.C

\*---------------------------------------------------------------------------*/

#ifndef writers_H
#define writers_H

#include "writer.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Only used internally
#define makeTypeSetWritersTypeName(typeWriter, dataType)                      \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(typeWriter< dataType >, 0)

// Sometimes used externally
#define makeSetWritersTypeName(typeWriter)                                    \
                                                                              \
    makeTypeSetWritersTypeName(typeWriter, scalar);                           \
    makeTypeSetWritersTypeName(typeWriter, vector);                           \
    makeTypeSetWritersTypeName(typeWriter, sphericalTensor);                  \
    makeTypeSetWritersTypeName(typeWriter, symmTensor);                       \
    makeTypeSetWritersTypeName(typeWriter, tensor);                           \
    makeTypeSetWritersTypeName(typeWriter, symmTensor4thOrder)

// Define type info for single dataType template instantiation (eg, vector)
#define makeSetWriterType(typeWriter, dataType)                               \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(typeWriter< dataType >, 0);           \
    addTemplatedToRunTimeSelectionTable                                       \
    (                                                                         \
        writer, typeWriter, dataType, word                                    \
    )


// Define type info for scalar, vector etc. instantiations
#define makeSetWriters(typeWriter)                                            \
                                                                              \
    makeSetWriterType(typeWriter, scalar);                                    \
    makeSetWriterType(typeWriter, vector);                                    \
    makeSetWriterType(typeWriter, sphericalTensor);                           \
    makeSetWriterType(typeWriter, symmTensor);                                \
    makeSetWriterType(typeWriter, tensor);                                    \
    makeSetWriterType(typeWriter, symmTensor4thOrder)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
