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
    tetPolyPatchMapperCellDecomp

Description
    Mapping class for a tetPolyPatchFieldCellDecomp

SourceFiles
    tetPolyPatchMapperCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyPatchMapperCellDecomp_H
#define tetPolyPatchMapperCellDecomp_H

#include "PointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyPatchCellDecomp;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                Class tetPolyPatchMapperCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyPatchMapperCellDecomp
:
    public PointPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const tetPolyPatchCellDecomp& patch_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;


    // Demand-driven private data

        //- Is the mapping direct
        mutable bool* directPtr_;

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyPatchMapperCellDecomp
        (
            const tetPolyPatchMapperCellDecomp&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPolyPatchMapperCellDecomp&);


        //- Calculate addressing for mapping with inserted cells
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Static data members

    // Constructors

        //- Construct from mappers
        tetPolyPatchMapperCellDecomp
        (
            const tetPolyPatchCellDecomp& patch,
            const mapPolyMesh& meshMap
        );


    // Destructor

        virtual ~tetPolyPatchMapperCellDecomp();


    // Member Functions

        //- Return size
        virtual label size() const;

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const
        {
            return true;
        }

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
