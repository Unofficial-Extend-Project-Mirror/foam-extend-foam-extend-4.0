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
    tetPolyPatchMapperFaceDecomp

Description
    Mapping class for a tetPolyPatchFieldFaceDecomp

SourceFiles
    tetPolyPatchMapperFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatchMapperFaceDecomp_H
#define tetPolyPatchMapperFaceDecomp_H

#include "PointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyPatchFaceDecomp;
class mapPolyMesh;
class pointMapper;
class faceMapper;

/*---------------------------------------------------------------------------*\
                Class tetPolyPatchMapperFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyPatchMapperFaceDecomp
:
    public PointPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const tetPolyPatchFaceDecomp& patch_;

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
        tetPolyPatchMapperFaceDecomp
        (
            const tetPolyPatchMapperFaceDecomp&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPolyPatchMapperFaceDecomp&);


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Static data members

    // Constructors

        //- Construct from mappers
        tetPolyPatchMapperFaceDecomp
        (
            const tetPolyPatchFaceDecomp& patch,
            const mapPolyMesh& meshMap,
            const pointMapper& pMapper,
            const faceMapper& fMapper
        );


    // Destructor

        virtual ~tetPolyPatchMapperFaceDecomp();


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
