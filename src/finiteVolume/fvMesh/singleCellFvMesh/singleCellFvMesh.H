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
    Foam::singleCellFvMesh

Description
    fvMesh as subset of other mesh. Consists of one cell and all original
    bounday faces. Useful when manipulating boundary data. Single internal
    cell only needed to be able to manipulate in a standard way.

SourceFiles
    singleCellFvMesh.C
    singleCellFvMeshInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef singleCellFvMesh_H
#define singleCellFvMesh_H

#include "fvPatchFieldMapper.H"
#include "fvMesh.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class singleCellFvMesh Declaration
\*---------------------------------------------------------------------------*/

class singleCellFvMesh
:
    public fvMesh
{
    // Private data

        const labelListIOList patchFaceAgglomeration_;

        //- From patch faces back to agglomeration or fine mesh
        labelListIOList patchFaceMap_;

        //- From fine mesh faces to coarse mesh
        labelIOList reverseFaceMap_;

        //- From coarse points back to original mesh
        labelIOList pointMap_;

        //- From fine points to coarse mesh
        labelIOList reversePointMap_;


    // Private Member Functions

        //- Calculate agglomerated mesh
        void agglomerateMesh(const fvMesh&, const labelListList&);


        //- Disallow default bitwise copy construct
        singleCellFvMesh(const singleCellFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const singleCellFvMesh&);


public:

        //- Patch field mapper class for non-agglomerated meshes
        class directPatchFieldMapper
        :
            public fvPatchFieldMapper
        {
            // Private data

                //- Mapping
                const UList<label>& directAddressing_;

                //- Size before mapping
                const label sizeBeforeMapping_;


        public:

                //- Construct given addressing
                directPatchFieldMapper
                (
                    const UList<label>& directAddressing,
                    const label sizeBeforeMapping
                )
                :
                    directAddressing_(directAddressing),
                    sizeBeforeMapping_(sizeBeforeMapping)
                {}

                virtual label size() const
                {
                    return directAddressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                virtual bool direct() const
                {
                    return true;
                }

                virtual const UList<label>& directAddressing() const
                {
                    return directAddressing_;
                }
        };

        //- Patch field mapper class for agglomerated meshes
        class agglomPatchFieldMapper
        :
            public fvPatchFieldMapper
        {
            // Private data

                //- Addressing
                const labelListList& addressing_;

                //- Weights
                const scalarListList& weights_;

                //- Size before mapping
                label sizeBeforeMapping_;


        public:

                //- Construct given addressing
                agglomPatchFieldMapper
                (
                    const labelListList& addressing,
                    const scalarListList& weights,
                    const label sizeBeforeMapping
                )
                :
                    addressing_(addressing),
                    weights_(weights),
                    sizeBeforeMapping_(sizeBeforeMapping)
                {}

                virtual label size() const
                {
                    return addressing_.size();
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                virtual bool direct() const
                {
                    return false;
                }

                virtual const labelListList& addressing() const
                {
                    return addressing_;
                }

                virtual const scalarListList& weights() const
                {
                    return weights_;
                }
        };



    // Constructors

        //- Construct from fvMesh and no agglomeration
        singleCellFvMesh(const IOobject& io, const fvMesh&);

        //- Construct from fvMesh and agglomeration of boundary faces.
        //  agglomeration is per patch, per patch face index the agglomeration
        //  the face goes into.
        singleCellFvMesh
        (
            const IOobject& io,
            const fvMesh&,
            const labelListList& patchFaceAgglomeration
        );

        //- Read from IOobject
        singleCellFvMesh(const IOobject& io);

    // Member Functions

        bool agglomerate() const
        {
            return patchFaceAgglomeration_.size() > 0;
        }

        //- From patchFace on this back to original mesh or agglomeration
        const labelListList& patchFaceMap() const
        {
            return patchFaceMap_;
        }

        //- From point on this back to original mesh
        const labelList& pointMap() const
        {
            return pointMap_;
        }

        //- From face on original mesh to face on this
        const labelList& reverseFaceMap() const
        {
            return reverseFaceMap_;
        }

        //- From point on original mesh to point on this (or -1 for removed
        //  points)
        const labelList& reversePointMap() const
        {
            return reversePointMap_;
        }

        //- Map volField. Internal field set to average, patch fields straight
        //  copies.
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh> >
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "singleCellFvMeshInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
