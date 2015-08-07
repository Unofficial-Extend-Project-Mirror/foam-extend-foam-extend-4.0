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
    tetPointFieldReconstructor

Description
    Tet point field reconstructor.

SourceFiles
    tetPointFieldReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef tetPointFieldReconstructor_H
#define tetPointFieldReconstructor_H

#include "PtrList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "tetPolyMesh.H"
#include "tetPointMesh.H"
#include "elementMesh.H"
#include "IOobjectList.H"
#include "tetPolyPatchFields.H"
#include "elementPatchFields.H"
#include "PointPatchFieldMapper.H"
#include "GeometricField.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class tetPointFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class tetPointFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        const tetPolyMesh& mesh_;

        //- List of processor meshes
        const PtrList<tetPolyMesh>& procMeshes_;

        //- List of processor point addressing lists
        const PtrList<labelIOList>& pointProcAddressing_;

        //- List of processor face addressing lists
        const PtrList<labelIOList>& faceProcAddressing_;

        //- List of processor cell addressing lists
        const PtrList<labelIOList>& cellProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPointFieldReconstructor
        (
            const tetPointFieldReconstructor&
        );

        //- Disallow default bitwise assignment
        void operator=(const tetPointFieldReconstructor&);

        //- Processor-to-global addressing for all points
        labelList procAddressing(const label procNo) const;

        //- Patch-to-patch addressing
        labelList procPatchAddressing
        (
            const labelList& procToGlobalAddr,
            const label procNo,
            const label patchNo
        ) const;


public:

        class tetPolyPatchFieldReconstructor
        :
            public PointPatchFieldMapper
        {
            label size_;
            label sizeBeforeMapping_;

        public:

            // Constructors

                //- Construct given addressing
                tetPolyPatchFieldReconstructor
                (
                    const label size,
                    const label sizeBeforeMapping
                )
                :
                    size_(size),
                    sizeBeforeMapping_(sizeBeforeMapping)
                {}


            // Member functions

                virtual label size() const
                {
                    return size_;
                }

                virtual label sizeBeforeMapping() const
                {
                    return sizeBeforeMapping_;
                }

                virtual bool direct() const
                {
                    return true;
                }

                virtual const unallocLabelList& directAddressing() const
                {
                    return unallocLabelList::null();
                }
        };


    // Constructors

        //- Construct from components
        tetPointFieldReconstructor
        (
            tetPolyMesh& mesh,
            const PtrList<tetPolyMesh>& procMeshes,
            const PtrList<labelIOList>& pointProcAddressing,
            const PtrList<labelIOList>& faceProcAddressing,
            const PtrList<labelIOList>& cellProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Reconstruct point field
        template<class Type>
        tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
        reconstructTetPointField
        (
            const IOobject& fieldIoObject
        );

        //- Reconstruct element field
        template<class Type>
        tmp<GeometricField<Type, elementPatchField, elementMesh> >
        reconstructElementField
        (
            const IOobject& fieldIoObject
        );

        //- Reconstruct and write all point fields
        template<class Type>
        void reconstructTetPointFields
        (
            const IOobjectList& objects
        );

        //- Reconstruct and write all element fields
        template<class Type>
        void reconstructElementFields
        (
            const IOobjectList& objects
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "tetPointFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
