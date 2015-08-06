/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    processorBoundaryPatch

Description
    Writes information for a given processor patch

SourceFiles
    processorBoundaryPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorBoundaryPatch_H
#define processorBoundaryPatch_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryPatchBase.H"
#include "typeInfo.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class processorBoundaryPatch
:   public boundaryPatchBase
{
    // private data

        label myProcNo_;

        label neighbProcNo_;

public:

    // Runtime type information
        TypeName("processor");

    // Construct from components
        processorBoundaryPatch
        (
            const word& name,
            const word& type,
            const label nFaces,
            const label startFace,
            const label myProcNo,
            const label neighbProcNo
        );

        processorBoundaryPatch(const word& name, const dictionary&);

    // Construct and return a clone
        virtual autoPtr<boundaryPatchBase> clone() const
        {
            return autoPtr<boundaryPatchBase>
            (
                new processorBoundaryPatch
                (
                    patchName(),
                    patchType(),
                    patchSize(),
                    patchStart(),
                    myProcNo_,
                    neighbProcNo_
                )
            );
        }

        virtual autoPtr<boundaryPatchBase> clone
        (
            const processorBoundaryPatch& pp
        ) const
        {
            return autoPtr<boundaryPatchBase>
            (
                new processorBoundaryPatch
                (
                    pp.patchName(),
                    pp.patchType(),
                    pp.patchSize(),
                    pp.patchStart(),
                    pp.myProcNo_,
                    pp.neighbProcNo_
                )
            );
        }

    // Member functions
        //- return the owner processor
        inline label myProcNo() const
        {
            return myProcNo_;
        }

        //- return the neighbour processor
        inline label neiProcNo() const
        {
            return neighbProcNo_;
        }

        //- check if the processor is the owner of the interface
        inline bool owner() const
        {
            if( myProcNo_ < neighbProcNo_ )
                return true;

            return false;
        }

        //- Return as dictionary of entries
        dictionary dict() const;

    // Write
        //- write to Ostream
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&) const;

    // Member operators
        Ostream& operator<<(Ostream&) const;

        Istream& operator>>(Istream&);

        void operator=(const processorBoundaryPatch&);

        bool operator!=(const processorBoundaryPatch&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
