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
    topoPatchMapper

Description
    This object provides mapping and fill-in information for patch data
    between the two meshes after the topological change.  It is
    constructed from mapPolyMesh.

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    topoPatchMapper.C

\*---------------------------------------------------------------------------*/

#ifndef topoPatchMapper_H
#define topoPatchMapper_H

#include "morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class topoMapper;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                    Class topoPatchMapper Declaration
\*---------------------------------------------------------------------------*/

class topoPatchMapper
:
    public fvPatchFieldMapper
{
    // Private data

        //- Reference to patch
        const fvPatch& patch_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Reference to the topoMapper
        const topoMapper& tMapper_;

        //- Is the mapping direct
        bool direct_;

    // Demand-driven private data

        //- Direct addressing
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing
        mutable labelListList* interpolationAddrPtr_;

        //- Inverse-distance weights
        mutable scalarListList* weightsPtr_;

        //- Inserted faces
        mutable labelList* insertedFaceLabelsPtr_;

        //- Inserted face index map
        mutable labelList* insertedFaceIndexMapPtr_;

        //- Inserted face addressing
        mutable labelListList* insertedFaceAddressingPtr_;

        //- Interpolation areas
        mutable List<scalarField>* areasPtr_;

        //- Interpolation centres
        mutable List<vectorField>* centresPtr_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        topoPatchMapper(const topoPatchMapper&);

        //- Disallow default bitwise assignment
        void operator=(const topoPatchMapper&);

        //- Calculate the insertedFace addressing list
        void calcInsertedFaceAddressing() const;

        //- Calculate addressing for mapping
        void calcAddressing() const;

        //- Calculate inverse-distance weights for interpolative mapping
        void calcInverseDistanceWeights() const;

        //- Calculate intersection weights for conservative mapping
        void calcIntersectionWeightsAndCentres() const;

        //- Return intersection area weights
        const List<scalarField>& intersectionWeights() const;

        //- Return intersection area centres
        const List<vectorField>& intersectionCentres() const;

        //- Clear out local storage
        void clearOut();

public:

    // Constructors

        //- Construct from mapPolyMesh
        topoPatchMapper
        (
            const fvPatch& patch,
            const mapPolyMesh& mpm,
            const topoMapper& mapper
        );

    // Destructor

        virtual ~topoPatchMapper();

    // Member Functions

        //- Return the polyPatch size
        label patchSize() const;

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

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const;

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;

        //- Return addressing map for inserted faces
        const labelList& insertedObjectMap() const;

        //- Return addressing for inserted faces
        const labelListList& insertedFaceAddressing() const;

        //- Map the patch field
        template <class Type>
        void mapPatchField
        (
            const word& fieldName,
            Field<Type>& pF
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#include "topoPatchMapper.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
