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
    boundaryPatchBase

Description
    Writes information for a given boundary patch

SourceFiles
    boundaryPatchBase.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryPatchBase_H
#define boundaryPatchBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "word.H"
#include "label.H"
#include "objectRegistry.H"
#include "typeInfo.H"
#include "dictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class boundaryPatchBase
{
protected:
    
    // Protected data
        word name_;
    
        word type_;
    
        label nFaces_;
    
        label startFace_;

public:
    
    //- Runtime type information
    TypeName("boundaryPatchBase");

    // Constructors
    
        boundaryPatchBase
        (
            const word&,
            const word&,
            const label,
            const label
        );

        boundaryPatchBase(const word& name, const dictionary&);

    // Destructor

        virtual ~boundaryPatchBase()
        {}

    // Declare run-time constructor selection table
            
        declareRunTimeSelectionTable
        (
            autoPtr,
            boundaryPatchBase,
            dictionary,
            (
                const word& name,
                const dictionary& dict
            ),
            (name, dict)
        );

    // Selectors
        
        static autoPtr<boundaryPatchBase> New
        (
            const word& name,
            const dictionary& dict
        );
        
        static autoPtr<boundaryPatchBase> New
        (
            Istream&
        );
        
    // Construct and return a clone
    
        autoPtr<boundaryPatchBase> clone() const
        {
            return New
            (
                name_,
                this->dict()
            );
        }
        
    //- Write
        virtual void write(Ostream&) const = 0;

        //- Write dictionary
        virtual void writeDict(Ostream&) const = 0;

    // Member functions

        inline const word& patchName() const
        {
            return name_;
        }
        
        inline word& patchName()
        {
            return name_;
        }
        
        inline const word& patchType() const
        {
            return type_;
        }
        
        inline word& patchType()
        {
            return type_;
        }
        
        inline label patchStart() const
        {
            return startFace_;
        }
        
        inline label& patchStart()
        {
            return startFace_;
        }
        
        inline label patchSize() const
        {
            return nFaces_;
        }
        
        inline label& patchSize()
        {
            return nFaces_;
        }
        
        //- Return as dictionary of entries
        virtual dictionary dict() const = 0;
        
    // Member operators
        friend Ostream& operator<<(Ostream&, const boundaryPatchBase&);
        virtual Istream& operator>>(Istream&) = 0;
        //virtual void operator=(const boundaryPatchBase&) = 0;
        //virtual bool operator!=(const boundaryPatchBase&) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
