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
        virtual void calcTransforms();

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


    // Destructor

        virtual ~cyclicGgiPolyPatch();


    // Member functions

        //- Return shadow patch
        const cyclicGgiPolyPatch& cyclicShadow() const;

        //- Return separation offset
        const vector& separationOffset() const
        {
            return separationOffset_;
        }

        //- Return rotation axis
        const vector& rotationAxis() const
        {
            return rotationAxis_;
        }

        //- Return rotation angle
        const scalar& rotationAngle() const
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
