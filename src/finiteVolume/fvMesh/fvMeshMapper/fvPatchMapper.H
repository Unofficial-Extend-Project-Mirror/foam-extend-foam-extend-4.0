/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::fvPatchMapper

Description
    Mapping class for a fvPatchField.

SourceFiles
    fvPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef fvPatchMapper_H
#define fvPatchMapper_H

#include "fvPatchFieldMapper.H"
#include "faceMapper.H"
#include "fvPatch.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvPatch;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class fvPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class fvPatchMapper
:
    public fvPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to face mapper
        const faceMapper& faceMap_;

        //- Size before mapping
        const label sizeBeforeMapping_;


    // Demand-driven private data

        //- Is the mapping direct
        mutable bool* directPtr_;

        //- Direct addressing (only one form of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one form of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvPatchMapper
        (
            const fvPatchMapper&
        );

        //- Disallow default bitwise assignment
        void operator=(const fvPatchMapper&);


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from mappers
        fvPatchMapper
        (
            const fvPatch& patch,
            const faceMapper& faceMap
        );


    // Destructor

        virtual ~fvPatchMapper();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return patch_.size();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return faceMap_.direct();
        }

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
