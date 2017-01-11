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
    cyclicGgiPolyPatch

Description
    Non-matching cyclic generalised grid interface (GGI) patch.

Author
    Martin Beaudoin, Hydro-Quebec, (2008)

Contributor:
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    cyclicGgiPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicGgiPolyPatch_H
#define cyclicGgiPolyPatch_H

#include "ggiPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cyclicGgiPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class cyclicGgiPolyPatch
:
    public ggiPolyPatch
{
    // Private data

        // Translation/rotation parameter for a cyclic

            //- Separation offset
            vector separationOffset_;

            //- Rotation axis
            vector rotationAxis_;

            //- Rotation angle
            scalar rotationAngle_;


    // Private member functions

        //- Calculate cyclic transforms (rotation and translation)
        //  Virtual over-ride for base GGI patch.  HJ, 14/Jan/2009
        virtual void calcTransforms() const;

        //- Check definition: angles and offsets
        void checkDefinition() const;


public:

    //- Runtime type information
    TypeName("cyclicGgi");


    // Constructors

        //- Construct from components
        cyclicGgiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from components
        cyclicGgiPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& shadowName,
            const word& zoneName,
            const bool bridgeOverlap,
            const vector& separationOffset,
            const vector& rotationAxis,
            const scalar rotationAngle
        );

        //- Construct from dictionary
        cyclicGgiPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        cyclicGgiPolyPatch
        (
            const cyclicGgiPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        cyclicGgiPolyPatch
        (
            const cyclicGgiPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new cyclicGgiPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new cyclicGgiPolyPatch
                (
                    refCast<const cyclicGgiPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~cyclicGgiPolyPatch();


    // Member functions

        //- Return shadow patch
        const cyclicGgiPolyPatch& cyclicShadow() const;

        //- Return separation offset
        const vector& separationOffset() const
        {
            return separationOffset_;
        }

        //- Return access to separation offset
        vector& separationOffset()
        {
            return separationOffset_;
        }

        //- Return rotation axis
        const vector& rotationAxis() const
        {
            return rotationAxis_;
        }

        //- Return access tp rotation axis
        vector& rotationAxis()
        {
            return rotationAxis_;
        }

        //- Return rotation angle
        const scalar& rotationAngle() const
        {
            return rotationAngle_;
        }

        //- Return access to rotation angle
        scalar& rotationAngle()
        {
            return rotationAngle_;
        }

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
