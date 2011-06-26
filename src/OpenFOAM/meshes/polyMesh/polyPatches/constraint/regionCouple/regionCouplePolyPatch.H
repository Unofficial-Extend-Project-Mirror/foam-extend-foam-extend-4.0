/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-7 H. Jasak All rights reserved
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
    regionCouplePolyPatch

Description
    Region couple patch.  Used for multi-region conjugate simulations

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    regionCouplePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef regionCouplePolyPatch_H
#define regionCouplePolyPatch_H

#include "coupledPolyPatch.H"
#include "patchToPatchInterpolation.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                     Class regionCouplePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCouplePolyPatch
:
    public coupledPolyPatch
{
    // Private data

        //- Shadow region name
        const word shadowRegionName_;

        //- Shadow patch name
        const word shadowPatchName_;

        //- Are the regions attached
        mutable Switch attached_;

        //- Are the region attached walls
        Switch isWall_;

        //- Shadow patch index.  Delayed evaluation for construction
        mutable label shadowIndex_;


        //- Patch-to-patch interpolation
        mutable patchToPatchInterpolation* patchToPatchPtr_;

        //- Reconstructed patch neighbour cell centres
        mutable vectorField* reconFaceCellCentresPtr_;


    // Private member functions

        //- Calculate interpolation
        void calcInterpolation() const;


protected:

    // Protected Member functions

        //- Initialise the calculation of the patch addressing
        virtual void initAddressing();

        //- Calculate the patch addressing
        virtual void calcAddressing();

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


        //- Return non-constant access to shadow patch
        regionCouplePolyPatch& shadow();

        //- Return reference to patch-to-patch interpolation
        const patchToPatchInterpolation& patchToPatch() const;

        //- Calculate reconstructed cell centres
        void calcReconFaceCellCentres() const;

        //- Clear geometry
        void clearGeom() const;

        //- Clear out
        void clearOut() const;


public:

    //- Runtime type information
    TypeName("regionCouple");


    // Constructors

        //- Construct from components
        regionCouplePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& shadowRegionName,
            const word& shadowPatchName,
            const bool attached,
            const bool attachedWall
        );

        //- Construct from dictionary
        regionCouplePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        regionCouplePolyPatch
        (
            const regionCouplePolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        regionCouplePolyPatch
        (
            const regionCouplePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new regionCouplePolyPatch(*this, bm));
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
                new regionCouplePolyPatch
                (
                    refCast<const regionCouplePolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~regionCouplePolyPatch();


    // Member functions

        //- Return true if coupled
        virtual bool coupled() const;

        //- Return shadow region name
        const word& shadowRegionName() const
        {
            return shadowRegionName_;
        }

        //- Return shadow patch name
        const word& shadowPatchName() const
        {
            return shadowPatchName_;
        }

        //- Return shadow region
        const polyMesh& shadowRegion() const;

        //- Return shadow patch index
        int shadowIndex() const;

        //- Return true if regions are attached
        bool attached() const
        {
            return attached_;
        }

        // Return true if patch is considered a wall
        bool isWall() const
        {
            return isWall_;
        }

        //- Attach regions
        void attach() const;

        //- Attach regions
        void detach() const;

        //- Return shadow patch
        const regionCouplePolyPatch& shadow() const;

        //- Interpolate face field
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& pf) const
        {
            return patchToPatch().faceInterpolate(pf);
        }

        template<class Type>
        tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const
        {
            return patchToPatch().faceInterpolate(tpf);
        }

        //- Return reconstructed cell centres
        const vectorField& reconFaceCellCentres() const;

        //- Force calculation of transformation tensors
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr
        ) const;

        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(const primitivePatch&) const;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation: for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        //- Synchronise communications of ordering for primitivePatch
        //  Used in cases when no topological change happens locally,
        //  but is happening on other processors
        virtual void syncOrder() const;


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
