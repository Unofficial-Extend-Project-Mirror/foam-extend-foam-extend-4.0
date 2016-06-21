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
    Foam::pointPatchMapper

Description
    Mapping class for a pointPatchField.

SourceFiles
    pointPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatchMapper_H
#define pointPatchMapper_H

#include "pointMapper.H"
#include "PointPatchFieldMapper.H"
#include "pointPatch.H"
#include "primitiveFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointPatch;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class pointPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class pointPatchMapper
:
    public PointPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const pointPatch& patch_;

        //- Reference to point mapper for internal field
        const morphFieldMapper& pointMapper_;

        //- Reference to mapping data
        const mapPolyMesh& mpm_;

        //- Size before mapping
        const label sizeBeforeMapping_;


    // Demand-driven private data

        //- Direct addressing (only one form of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one form of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointPatchMapper
        (
            const pointPatchMapper&
        );

        //- Disallow default bitwise assignment
        void operator=(const pointPatchMapper&);


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Static data members

    // Constructors

        //- Construct from mappers
        pointPatchMapper
        (
            const pointPatch& patch,
            const pointMapper& pointMap,
            const mapPolyMesh& mpm
        );


    // Destructor

        virtual ~pointPatchMapper();


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
            return pointMapper_.direct();
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
