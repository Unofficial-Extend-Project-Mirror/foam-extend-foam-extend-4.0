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
    ggiFvPatch

Description
    Generalized grid interface patch.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    ggiFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef ggiFvPatch_H
#define ggiFvPatch_H

#include "coupledFvPatch.H"
#include "ggiLduInterface.H"
#include "ggiPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ggiFvPatch Declaration
\*---------------------------------------------------------------------------*/

class ggiFvPatch
:
    public coupledFvPatch,
    public ggiLduInterface
{
    // Private Data

        //- Reference to polyPatch
        const ggiPolyPatch& ggiPolyPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const;

        //- Make patch face non-orthogonality correction vectors
        virtual void makeCorrVecs(vectorField&) const;


public:

    //- Runtime type information
    TypeName(ggiPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        ggiFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            ggiPolyPatch_(refCast<const ggiPolyPatch>(patch))
        {}


    //- Destructor
    virtual ~ggiFvPatch();


    // Member functions

        // Communication support

            //- Return communicator used for parallel communication
            virtual int comm() const
            {
                return ggiPolyPatch_.comm();
            }

            //- Return message tag used for sending
            virtual int tag() const
            {
                return ggiPolyPatch_.tag();
            }


        // Access

            //- Return shadow patch
            const ggiFvPatch& shadow() const;

            //- Use bridging to fix overlap error in interpolation
            bool bridgeOverlap() const
            {
                return ggiPolyPatch_.bridgeOverlap();
            }

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;


        // Interpolation

            //- Interpolate face field
            template<class Type>
            tmp<Field<Type> > interpolate(const Field<Type>& pf) const
            {
                return ggiPolyPatch_.interpolate(pf);
            }

            template<class Type>
            tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const
            {
                return ggiPolyPatch_.interpolate(tpf);
            }

            //- Bridge interpolated face field for uncovered faces
            template<class Type>
            void bridge
            (
                const Field<Type>& bridgeField,
                Field<Type>& ff
            ) const
            {
                return ggiPolyPatch_.bridge(bridgeField, ff);
            }


        // Interface transfer functions

            //- Is this the master side?
            virtual bool master() const;

            //- Is this the fine level?
            virtual bool fineLevel() const;

            //- Return shadow patch index
            virtual label shadowIndex() const;

            //- Return shadow interface
            virtual const ggiLduInterface& shadowInterface() const;

            //- Return interface size
            virtual label interfaceSize() const;

            //- Return zone size
            virtual label zoneSize() const;

            //- Return zone addressing
            virtual const labelList& zoneAddressing() const;

            //- Return addressing.  Master side returns own addressing and
            //  slave side returns addressing from master
            virtual const labelListList& addressing() const;

            //- Is the patch localised on a single processor
            virtual bool localParallel() const;

            //- Return mapDistribute
            virtual const mapDistribute& map() const;

            //- Return weights.  Master side returns own weights and
            //  slave side returns weights from master
            virtual const scalarListList& weights() const;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return coupledFvPatch::forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return coupledFvPatch::reverseT();
            }

            //- Expand addressing to zone
            //  Used in optimised AMG coarsening
            virtual void expandAddrToZone(labelField&) const;


            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Initialise transfer of internal field adjacent to the interface
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& internalData
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
