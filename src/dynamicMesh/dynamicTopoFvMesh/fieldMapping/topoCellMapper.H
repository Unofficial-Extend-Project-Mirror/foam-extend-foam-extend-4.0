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
    topoCellMapper

Description
    This object provides mapping and fill-in information for internal
    cell data between the two meshes after topological changes.
    It is constructed from mapPolyMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    topoCellMapper.C

\*---------------------------------------------------------------------------*/

#ifndef topoCellMapper_H
#define topoCellMapper_H

#include "topoMapper.H"
#include "morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class topoCellMapper Declaration
\*---------------------------------------------------------------------------*/

class topoCellMapper
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

        //- Inserted cells
        mutable labelList* insertedCellLabelsPtr_;

        //- Interpolation volumes
        mutable List<scalarField>* volumesPtr_;

        //- Interpolation centres
        mutable List<vectorField>* centresPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoCellMapper(const topoCellMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoCellMapper&);

        //- Calculate addressing for interpolative mapping
        void calcAddressing() const;

        //- Calculate inverse-distance weights for interpolative mapping
        void calcInverseDistanceWeights() const;

        //- Calculate intersection weights for conservative mapping
        void calcIntersectionWeightsAndCentres() const;

        //- Return intersection volume weights
        const List<scalarField>& intersectionWeights() const;

        //- Return intersection volume centres
        const List<vectorField>& intersectionCentres() const;

        //- Clear out local storage
        void clearOut();

public:

    // Constructors

        //- Construct from mapPolyMesh
        topoCellMapper
        (
            const mapPolyMesh& mpm,
            const topoMapper& mapper
        );

    // Destructor

        virtual ~topoCellMapper();

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

        //- Are there any inserted cells
        virtual bool insertedObjects() const;

        //- Return list of inserted cells
        virtual const labelList& insertedObjectLabels() const;

        //- Conservatively map the internal field
        template <class Type, class gradType>
        void mapInternalField
        (
            const word& fieldName,
            const Field<gradType>& gF,
            Field<Type>& iF
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "topoCellMapperTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
