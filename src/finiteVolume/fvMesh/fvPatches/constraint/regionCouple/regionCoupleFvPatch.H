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
    regionCoupleFvPatch

Description
    Region couple patch coupled two mesh regions for a solution in the same
    linear system.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    regionCoupleFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef regionCoupleFvPatch_H
#define regionCoupleFvPatch_H

#include "coupledFvPatch.H"
#include "regionCoupleLduInterface.H"
#include "regionCouplePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                     Class regionCoupleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCoupleFvPatch
:
    public coupledFvPatch,
    public regionCoupleLduInterface
{
    // Private Data

        //- Reference to polyPatch
        const regionCouplePolyPatch& rcPolyPatch_;

        //- Transfer buffer
        mutable labelField transferBuffer_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const;

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const;

        //- Return contents of transfer buffer
        const labelField& transferBuffer() const
        {
            return transferBuffer_;
        }


public:

    //- Runtime type information
    TypeName(regionCouplePolyPatch::typeName_());


    // Constructors

        //- Construct from components
        regionCoupleFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            rcPolyPatch_(refCast<const regionCouplePolyPatch>(patch)),
            transferBuffer_()
        {}


    // Destructor

        virtual ~regionCoupleFvPatch()
        {}


    // Member functions

        //- Return true if coupled
        virtual bool coupled() const;

        //- Return shadow patch index
        int shadowIndex() const
        {
            return rcPolyPatch_.shadowIndex();
        }

        //- Return shadow region
        const fvMesh& shadowRegion() const;

        //- Return shadow patch
        const regionCoupleFvPatch& shadow() const;

        //- Interpolate face field
        template<class Type>
        tmp<Field<Type> > interpolate(const Field<Type>& pf) const
        {
            return rcPolyPatch_.interpolate(pf);
        }

        template<class Type>
        tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const
        {
            return rcPolyPatch_.interpolate(tpf);
        }

        //- Return delta (P to N) vectors across coupled patch
        virtual tmp<vectorField> delta() const;


        // Interface transfer functions

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
