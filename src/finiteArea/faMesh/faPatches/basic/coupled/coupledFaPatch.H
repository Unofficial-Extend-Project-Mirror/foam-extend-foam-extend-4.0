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
    coupledFaPatch

Author
    Zeljko Tukovic and Hrvoje Jasak

Description
    coupledFaPatch is an abstract base class for patches that couple regions
    of the computational domain e.g. cyclic, arbitrary interfaces, sliding
    interfaces and processor-processor links.

SourceFiles
    coupledFaPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFaPatch_H
#define coupledFaPatch_H

#include "lduInterface.H"
#include "faPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coupledFaPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFaPatch
:
    public lduInterface,
    public faPatch
{
    // Private data

        //- offset (distance) vector from one side of the couple to the other
        mutable vectorField separation_;

        //- Face transformation tensor
        mutable tensorField forwardT_;

        //- Neighbour-cell transformation tensor
        mutable tensorField reverseT_;


protected:

    // Protected Member Functions

        //- Make patch weighting factors
        virtual void makeWeights(scalarField&) const = 0;

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const = 0;

        //- Calculate the uniform transformation tensors
        void calcTransformTensors
        (
            const vector& Cf,
            const vector& Cr,
            const vector& nf,
            const vector& nr
        ) const;

        //- Calculate the transformation tensors
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr
        ) const;


public:

    //- Runtime type information
    TypeName("coupled");


    // Constructors

        //- Construct from components
        coupledFaPatch
        (
            const word& name,
            const labelList& edgeLabels,
            const label index,
            const faBoundaryMesh& bm,
            const label ngbPolyPatchIndex
        )
        :
            faPatch(name, edgeLabels, index, bm, ngbPolyPatchIndex)
        {}

        //- Construct from dictionary
        coupledFaPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const faBoundaryMesh& bm
        )
        :
            faPatch(name, dict, index, bm)
        {}


    // Destructor

        virtual ~coupledFaPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }


            //- Are the coupled planes separated
            bool separated() const
            {
                return separation_.size();
            }

            //- Return the offset (distance) vector from one side of the couple
            //  to the other
            const vectorField& separation() const
            {
                if (!separation_.size())
                {
                    FatalErrorIn("coupledFaPatch::separation() const")
                        << "Coupled patches are not separated"
                        << abort(FatalError);
                }

                return separation_;
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                if (!forwardT_.size())
                {
                    FatalErrorIn("coupledFaPatch::forwardT() const")
                        << "Coupled planes do not need transformation"
                        << abort(FatalError);
                }

                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                if (!reverseT_.size())
                {
                    FatalErrorIn("coupledFaPatch::forwardT() const")
                        << "Coupled planes do not need transformation"
                        << abort(FatalError);
                }

                return reverseT_;
            }

            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return forwardT_.size() == 0;
            }


        //- Initialise the calculation of the patch geometry
        virtual void initGeometry() = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry() = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(const pointField&) = 0;


        // Access functions for demand driven data

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const = 0;


        // Interface transfer functions

            //- Return faceCell addressing: lduInterface virtual function
            virtual const unallocLabelList& faceCells() const
            {
                return edgeFaces();
            }

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const = 0;

            //- Initialise interface data transfer
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const
            {}

            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const = 0;

            //- Initialise neighbour field transfer
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                unallocLabelList& iF
            ) const
            {}

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
