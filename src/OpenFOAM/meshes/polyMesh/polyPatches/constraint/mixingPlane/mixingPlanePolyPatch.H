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
    mixingPlanePolyPatch

Description
    MixingPlane patch.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlanePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlanePolyPatch_H
#define mixingPlanePolyPatch_H

#include "coupledPolyPatch.H"
#include "mixingPlaneInterpolation.H"
#include "transform.H"
#include "word.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mixingPlanePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class mixingPlanePolyPatch
:
    public coupledPolyPatch
{
    // Private data

        //- Shadow patch name
        const word shadowName_;

        //- Coordinate system for the mixing plane
        autoPtr<coordinateSystem> csPtr_;

        //- Type of mixing plane assembly algorithm
        MixingPlaneInterpolationName::assembly assemblyType_;

        //- Orientation of the mixing plane
        //  direction = axis in which ribbon width is measured
        //  span = axis in which ribbon length is measured
        MixingPlaneInterpolationName::orientation orientationType_;

        //- User-defined interpolation profile file name
        fileName userProfileFile_;


        // Demand-driven data

            //- Shadow patch index.  Delayed evaluation for construction
            mutable label shadowIndex_;

            //- Patch-to-patch interpolation
            mutable mixingPlaneInterpolation* patchToPatchPtr_;

            //- Reconstructed patch neighbour cell centres
            mutable vectorField* reconFaceCellCentresPtr_;


    // Private member functions

        //- Calculate interpolation
        virtual void calcPatchToPatch() const;

        //- Calculate reconstructed cell centres
        void calcReconFaceCellCentres() const;

        //- Force calculation of transformation tensors
        virtual void calcTransforms();

        //- Clear out
        void clearOut();


protected:

    // Protected Member functions

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


public:

    //- Runtime type information
    TypeName("mixingPlane");


    // Constructors

        //- Construct from components
        mixingPlanePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from components
        mixingPlanePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm,
            const word& shadowName,
            const coordinateSystem& cs,
            const mixingPlaneInterpolation::assembly assemblyType,
            const mixingPlaneInterpolation::orientation orientationType,
            const fileName& userProfileFile
        );

        //- Construct from dictionary
        mixingPlanePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );


        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new mixingPlanePolyPatch(*this, bm));
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
                new mixingPlanePolyPatch
                (
                    refCast<const mixingPlanePolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~mixingPlanePolyPatch();


    // Member functions

        //- Return shadow patch name
        const word& shadowName() const
        {
            return shadowName_;
        }

        //- Return shadow patch index
        label shadowIndex() const;

        //- Return shadow patch
        const mixingPlanePolyPatch& shadow() const;

        //- Is this the master side?
        bool master() const
        {
            return index() < shadowIndex();
        }

        //- Is this the slave side?
        bool slave() const
        {
            return !master();
        }

        //- Return coordinate system
        const coordinateSystem& cs() const;

        //- Return reference to patch-to-patch interpolation
        //  Used only for addressing
        const mixingPlaneInterpolation& patchToPatch() const;

        //- Interpolate face field: given field on a the shadow side,
        //  create an interpolated field on this side
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& pf) const;

        template<class Type>
        tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const;

        //- Interpolate face field from average: given field on a the
        //- master/slave side, create an interpolated field on the same side
        //- using averaged values
        template<class Type>
        tmp<Field<Type> > circumferentialAverage(const Field<Type>& pf) const;

        template<class Type>
        tmp<Field<Type> > circumferentialAverage(const tmp<Field<Type> >& tpf) const;

        //- Return reconstructed cell centres
        const vectorField& reconFaceCellCentres() const;

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

#ifdef NoRepository
#   include "mixingPlanePolyPatchTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
