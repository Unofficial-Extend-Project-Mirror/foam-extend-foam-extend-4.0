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
    Foam::pointFieldReconstructor

Description
    Point field reconstructor.

SourceFiles
    pointFieldReconstructor.C

\*---------------------------------------------------------------------------*/

#ifndef pointFieldReconstructor_H
#define pointFieldReconstructor_H

#include "pointMesh.H"
#include "pointFields.H"
#include "PointPatchFieldMapperPatchRef.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class pointFieldReconstructor Declaration
\*---------------------------------------------------------------------------*/

class pointFieldReconstructor
{
    // Private data

        //- Reconstructed mesh reference
        const pointMesh& mesh_;

        //- List of processor meshes
        const PtrList<pointMesh>& procMeshes_;

        //- List of processor point addressing lists
        const PtrList<labelIOList>& pointProcAddressing_;

        //- List of processor boundary addressing lists
        const PtrList<labelIOList>& boundaryProcAddressing_;

        //- Point patch addressing
        labelListListList patchPointAddressing_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointFieldReconstructor
        (
            const pointFieldReconstructor&
        );

        //- Disallow default bitwise assignment
        void operator=(const pointFieldReconstructor&);


public:

        class pointPatchFieldReconstructor
        :
            public PointPatchFieldMapper
        {
            label size_;
            label sizeBeforeMapping_;


        public:

            // Constructors

                //- Construct given size
                pointPatchFieldReconstructor
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
        pointFieldReconstructor
        (
            const pointMesh& mesh,
            const PtrList<pointMesh>& procMeshes,
            const PtrList<labelIOList>& pointProcAddressing,
            const PtrList<labelIOList>& boundaryProcAddressing
        );


    // Member Functions

        //- Reconstruct field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> >
        reconstructField(const IOobject& fieldIoObject);

        //- Reconstruct and write all fields
        template<class Type>
        void reconstructFields(const IOobjectList& objects);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pointFieldReconstructorReconstructFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
