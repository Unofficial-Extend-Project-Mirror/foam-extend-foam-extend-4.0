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
    overlapGgiFvPatch

Description
    Partial overlap generalised grid interface (GGI) patch.  Master side
    remains unchanged and the slave side will be copied as needed
    to pave the master patch surface.

    This implies that the master patch has got a larger angular pitch than
    the slave and that master and slave are aligned at one edge.
    Master and slave will specify the pitch, based on which the expansion
    of the master side will be performed and used in interpolation.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    overlapGgiFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGgiFvPatch_H
#define overlapGgiFvPatch_H

#include "ggiFvPatch.H"
#include "overlapGGILduInterface.H"
#include "overlapGgiPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class overlapGgiFvPatch Declaration
\*---------------------------------------------------------------------------*/

class overlapGgiFvPatch
:
    public coupledFvPatch,
    public overlapGGILduInterface
{
    // Private Data

        //- Reference to polyPatch
        const overlapGgiPolyPatch& overlapGgiPolyPatch_;


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
    TypeName(overlapGgiPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        overlapGgiFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            overlapGgiPolyPatch_(refCast<const overlapGgiPolyPatch>(patch))
        {}


    // Destructor

        virtual ~overlapGgiFvPatch()
        {}


    // Member functions

        // Access

            //- Is this the master side?
            bool master() const
            {
                return overlapGgiPolyPatch_.master();
            }

            //- Return shadow patch
            const overlapGgiFvPatch& shadow() const;

            //- Interpolate face field
            template<class Type>
            tmp<Field<Type> > interpolate(const Field<Type>& pf) const
            {
                return overlapGgiPolyPatch_.interpolate(pf);
            }

            template<class Type>
            tmp<Field<Type> > interpolate(const tmp<Field<Type> >& tpf) const
            {
                return overlapGgiPolyPatch_.interpolate(tpf);
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

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
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
