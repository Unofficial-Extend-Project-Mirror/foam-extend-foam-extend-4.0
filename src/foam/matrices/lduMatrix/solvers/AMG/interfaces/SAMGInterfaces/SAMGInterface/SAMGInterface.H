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
    Foam::SAMGInterface

Description
    Abstract base class for AMG selected interfaces.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    SAMGInterface.C
    newSAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef SAMGInterface_H
#define SAMGInterface_H

#include "lduInterface.H"
#include "autoPtr.H"
#include "lduPrimitiveMesh.H"
#include "crMatrix.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class SAMGInterface
:
    public lduInterface
{
    // Private data

        //- Reference to ldu addressing
        const lduPrimitiveMesh& lduMesh_;

        //- Reference to filtered (interface only) prolongation matrix
        const crMatrix& interfaceProlongation_;

        //- Reference to filtered neighbour prolongation matrix
        const crMatrix nbrInterfaceProlongation_;


protected:

    // Protected data

        //- Face-cell addressing.  Contains coarse level addressing
        labelField faceCells_;

        //- Fine addressing
        //  On SAMG interfaces, a single fine coefficient may contribute to
        //  multiple coarse coefficients using different weights.
        //  To handle this, a fine coefficient may be visited multiple times
        //  which is recorded in fineAddressing.
        //  For simple (matching) interfaces, fineAddressing_[i] = i
        //  HJ, 21/Jun/2016
        labelField fineAddressing_;

        //- Restrict addressing
        //  For each fine coefficient, list coarse cluster index it will be
        //  agglomerated into
        //  For cases where the fineAddressing is used, restrict addressing
        //  and weights are expanded to match multiple hits for a single
        //  fine coefficient, as dictated by fineAddressing
        //  HJ, 21/Jun/2016
        labelField restrictAddressing_;

        //- Fine level agglomeration weights
        scalarField restrictWeights_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        SAMGInterface(const SAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const SAMGInterface&);


public:

    //- Runtime type information
    TypeName("SAMGInterface");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            SAMGInterface,
            lduInterface,
            (
                const lduPrimitiveMesh& lduMesh,
                const crMatrix& interfaceProlongation,
                const lduInterfacePtrsList& coarseInterfaces,
                const lduInterface& fineInterface,
                const crMatrix& nbrInterfaceProlongation
            ),
            (
                lduMesh,
                interfaceProlongation,
                coarseInterfaces,
                fineInterface,
                nbrInterfaceProlongation
            )
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<SAMGInterface> New
        (
            const lduPrimitiveMesh& lduMesh,
            const crMatrix& interfaceProlongation,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const crMatrix& nbrInterfaceProlongation
        );


    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour select addressing
        SAMGInterface
        (
            const lduPrimitiveMesh& lduMesh,
            const crMatrix& interfaceProlongation,
            const crMatrix& nbrInterfaceProlongation
        )
        :
            lduMesh_(lduMesh),
            interfaceProlongation_(interfaceProlongation),
            nbrInterfaceProlongation_(nbrInterfaceProlongation)
        {}


    //- Destructor
    virtual ~SAMGInterface()
    {}


    // Member Functions

        // Access

            //- Return reference to addressing
            const lduMesh& ldu() const
            {
                return lduMesh_;
            }

            //- Return reference to filtered prolongation matrix
            const crMatrix& interfaceProlongation() const
            {
                return interfaceProlongation_;
            }

            //- Return reference to neighbour's filtered prolongation matrix
            const crMatrix& nbrInterfaceProlongation() const
            {
                return nbrInterfaceProlongation_;
            }


            //- Return local size
            virtual label size() const
            {
                return faceCells_.size();
            }

            //- Return interface size
            virtual label interfaceSize() const
            {
                return faceCells_.size();
            }

            //- Return local faceCell addressing
            virtual const unallocLabelList& faceCells() const
            {
                return faceCells_;
            }

            //- Return fine addressing: fine index for each coarse face
            const labelField& fineAddressing() const
            {
                return fineAddressing_;
            }

            //- Return restrict addressing
            const labelField& restrictAddressing() const
            {
                return restrictAddressing_;
            }

            //- Return fine level agglomeration weights
            const scalarField& restrictWeights() const
            {
                return restrictWeights_;
            }

            //- Return the interface internal field of the given field
            template<class Type>
            tmp<Field<Type> > interfaceInternalField
            (
                const UList<Type>& internalData
            ) const;

            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;

            //- Initialise transfer of prolongation matrix to the interface
            virtual void initProlongationTransfer
            (
                const Pstream::commsTypes commsType,
                const crMatrix& filteredP
            ) const = 0;

            //- Transfer and return prolongation matrix adjacent to
            //  the interface
            virtual autoPtr<crMatrix> prolongationTransfer
            (
                const Pstream::commsTypes commsType,
                const crMatrix& filteredP
            ) const = 0;


        // Coefficient selection

            //- Select the given fine-level coefficients and return
            virtual tmp<scalarField> selectCoeffs
            (
                const scalarField& fineCoeffs
            ) const;

            // Coefficient agglomeration functions for block coefficients
            // moved to BlockSAMGInterfaceField to allow
            // templating on virtual functions in block matrix agglomeration
            // HJ, 16/Mar/2016
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SAMGInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
