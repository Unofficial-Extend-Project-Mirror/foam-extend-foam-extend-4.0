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
    Foam::faAreaMapper

Description
    FA area mapper.

SourceFiles
    faAreaMapper.C

\*---------------------------------------------------------------------------*/

#ifndef faAreaMapper_H
#define faAreaMapper_H

#include "morphFieldMapper.H"
#include "faMesh.H"
#include "faceMapper.H"
#include "HashSet.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class faAreaMapper Declaration
\*---------------------------------------------------------------------------*/

class faAreaMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to mesh mapper
        const faMesh& mesh_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Are there any inserted (unmapped) faces
        bool insertedFaces_;

        //- Is the mapping direct
        bool direct_;


    // Demand-driven private data

        //- Old mesh size
        label sizeBeforeMapping_;

        //- New face labels after mapping
        mutable labelList* newFaceLabelsPtr_;

        //- New face labels after mapping
        mutable labelList* newFaceLabelsMapPtr_;


        //- Direct addressing (only one form of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one form of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted faces
        mutable labelList* insertedObjectLabelsPtr_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        faAreaMapper(const faAreaMapper&);

        //- Disallow default bitwise assignment
        void operator=(const faAreaMapper&);


        //- Calculate addressing
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from components
        faAreaMapper
        (
            const faMesh& mesh,
            const mapPolyMesh& mpm
        );


    // Destructor

        virtual ~faAreaMapper();


    // Member Functions

        //- Return new face labels
        const labelList& newFaceLabels() const;

        //- Return new face labels map
        //  For new faces return old face index if it exists
        //  If the face has been added, index will be -1
        const labelList& newFaceLabelsMap() const;

        //- Return size
        virtual label size() const
        {
            return newFaceLabels().size();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return sizeBeforeMapping_;
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return direct_;
        }

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const
        {
            return !insertedObjectLabels().empty();
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
