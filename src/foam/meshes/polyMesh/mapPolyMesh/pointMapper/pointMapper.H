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
    Foam::pointMapper

Description
    This object provides mapping and fill-in information for point data
    between the two meshes after the topological change.  It is
    constructed from mapPolyMesh.

SourceFiles
    pointMapper.C

\*---------------------------------------------------------------------------*/

#ifndef pointMapper_H
#define pointMapper_H

#include "morphFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapPolyMesh;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class pointMapper Declaration
\*---------------------------------------------------------------------------*/

class pointMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Reference to mapPolyMesh
        const mapPolyMesh& mpm_;

        //- Are there any inserted (unmapped) points
        bool insertedPoints_;

        //- Is the mapping direct
        bool direct_;


    // Demand-driven private data

        //- Direct addressing (only one form of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted points
        mutable labelList* insertedPointLabelsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pointMapper(const pointMapper&);

        //- Disallow default bitwise assignment
        void operator=(const pointMapper&);


        //- Calculate addressing for mapping with inserted points
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from mapPolyMesh
        pointMapper(const mapPolyMesh& mpm);


    // Destructor

        virtual ~pointMapper();


    // Member Functions

        //- Return size
        virtual label size() const;

        //- Return size before mapping
        virtual label sizeBeforeMapping() const;

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

        //- Are there any inserted points
        bool insertedObjects() const
        {
            return insertedPoints_;
        }

        //- Return list of inserted points
        const labelList& insertedObjectLabels() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
