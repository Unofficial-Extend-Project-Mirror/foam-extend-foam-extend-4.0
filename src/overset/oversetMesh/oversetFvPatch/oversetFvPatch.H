/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::oversetFvPatch

Description
    Overset patch.

SourceFiles
    oversetFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef oversetFvPatch_H
#define oversetFvPatch_H

#include "coupledFvPatch.H"
#include "oversetLduInterface.H"
#include "oversetPolyPatch.H"
#include "oversetFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class oversetMesh;

/*---------------------------------------------------------------------------*\
                        Class oversetFvPatch Declaration
\*---------------------------------------------------------------------------*/

class oversetFvPatch
:
    public coupledFvPatch,
    public oversetLduInterface
{
    // Private Data

        //- Reference to overset patch
        const oversetPolyPatch& oversetPatch_;


protected:

    // Protected Member functions

        //- Make patch weighting factors
        void makeWeights(fvsPatchScalarField&) const;

        //- Make patch face - neighbour cell distances
        void makeDeltaCoeffs(fvsPatchScalarField&) const;


public:

    //- Runtime type information
    TypeName(oversetPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        oversetFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            coupledFvPatch(patch, bm),
            oversetPatch_(refCast<const oversetPolyPatch>(patch))
        {}


    //- Destructor
    virtual ~oversetFvPatch()
    {}


    // Member functions

        // Access

            //- Return delta (P to N) vectors across coupled patch
            virtual tmp<vectorField> delta() const;

            //- Return access to overset mesh
            virtual const oversetMesh& overset() const;

            //- Const access to overset poly patch
            const oversetPolyPatch& oversetPatch() const;

            //- Return true: coupled patch
            virtual bool coupled() const
            {
                return true;
            }


        // Interface transfer functions

            //- Is this the master side?
            virtual bool master() const;

            //- Is this the fine level?
            virtual bool fineLevel() const;

            //- Return interface size. On the fine level, this is the size
            //  of internal field
            virtual label interfaceSize() const;

            //- Return acceptor cells
            virtual const labelList& acceptorCells() const;

            //- Return donor cells
            virtual const labelList& donorCells() const;

            //- Return donor cells processor (to which donor is sent)
            virtual const labelList& donorCellsProc() const;

            //- Return mapDistribute
            virtual const mapDistribute& map() const;

            //- Return the values of the given internal data adjacent to
            //  the interface as a field.  For overset, this is fringe cell data
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;

            //- Initialise interface data transfer (face field).  Dummy
            virtual void initTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Transfer and return neighbour field (face field).  Dummy
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Initialise neighbour field transfer
            // For overset, this is donor data
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& internalData
            ) const;

            //- Return neighbour field
            //  For overset, this is donor data
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
