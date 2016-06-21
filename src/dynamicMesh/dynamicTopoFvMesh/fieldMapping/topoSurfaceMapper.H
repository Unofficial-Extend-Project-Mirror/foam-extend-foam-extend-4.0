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
    topoSurfaceMapper

Description
    This object provides mapping and fill-in information for internal
    face data between the two meshes after topological changes.
    It is constructed from mapPolyMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    topoSurfaceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef topoSurfaceMapper_H
#define topoSurfaceMapper_H

#include "topoMapper.H"
#include "morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class topoSurfaceMapper Declaration
\*---------------------------------------------------------------------------*/

class topoSurfaceMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Reference to the topoMapper
        const topoMapper& tMapper_;

        //- Is the mapping direct
        bool direct_;

        //- Size before mapping
        mutable label sizeBeforeMapping_;

    // Demand-driven private data

        //- Direct addressing
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing
        mutable labelListList* interpolationAddrPtr_;

        //- Inverse-distance weights
        mutable scalarListList* weightsPtr_;

        //- Inserted faces
        mutable labelList* insertedFaceLabelsPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoSurfaceMapper(const topoSurfaceMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoSurfaceMapper&);

        //- Calculate the insertedFaceLabels list
        void calcInsertedFaceLabels() const;

        //- Calculate addressing for mapping
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();

public:

    // Constructors

        //- Construct from mapPolyMesh
        topoSurfaceMapper
        (
            const mapPolyMesh& mpm,
            const topoMapper& mapper
        );

    // Destructor

        virtual ~topoSurfaceMapper();

    // Member Functions

        //- Return size
        virtual label size() const;

        //- Return size before mapping
        virtual label sizeBeforeMapping() const;

        //- Is the mapping direct
        virtual bool direct() const;

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

        //- Return interpolation addressing
        virtual const labelListList& addressing() const;

        //- Return weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const;

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;

        //- Return flux flip map
        const labelHashSet& flipFaceFlux() const;

        //- Map the internal field
        template <class Type>
        void mapInternalField
        (
            const word& fieldName,
            Field<Type>& iF
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "topoSurfaceMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
