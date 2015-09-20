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
    tetPolyPatchMapper

Description
    Mapping class for a tetPolyPatchField

SourceFiles
    tetPolyPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatchMapper_H
#define tetPolyPatchMapper_H

#include "PointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyPatch;
class mapPolyMesh;
class pointMapper;
class faceMapper;

/*---------------------------------------------------------------------------*\
                Class tetPolyPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class tetPolyPatchMapper
:
    public PointPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const tetPolyPatch& patch_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Reference to point mapper
        const pointMapper& pMapper_;

        //- Reference to face mapper
        const faceMapper& fMapper_;


    // Demand-driven private data

        //- Is the mapping direct
        mutable bool* directPtr_;

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyPatchMapper
        (
            const tetPolyPatchMapper&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPolyPatchMapper&);


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Static data members

    // Constructors

        //- Construct from mappers
        tetPolyPatchMapper
        (
            const tetPolyPatch& patch,
            const mapPolyMesh& meshMap,
            const pointMapper& pMapper,
            const faceMapper& fMapper
        );


    // Destructor

        virtual ~tetPolyPatchMapper();


    // Member Functions

        //- Return size
        virtual label size() const;

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const;

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolation weights
        virtual const scalarListList& weights() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
