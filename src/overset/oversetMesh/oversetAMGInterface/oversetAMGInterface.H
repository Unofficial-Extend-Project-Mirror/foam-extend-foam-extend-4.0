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
    Foam::oversetAMGInterface

Description
    AMG agglomerated overset interface.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    oversetAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef oversetAMGInterface_H
#define oversetAMGInterface_H

#include "AMGInterface.H"
#include "oversetLduInterface.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class oversetAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class oversetAMGInterface
:
    public AMGInterface,
    virtual public oversetLduInterface
{
    // Private data

        //- Reference tor the oversetLduInterface from which this is
        //  agglomerated
        const oversetLduInterface& fineOversetInterface_;

        //- Coarse interface size
        label interfaceSize_;

        // Local acceptor data

            // Local acceptor index is stored in faceCells_

            //- Donor cell index (can be on different processor)
            //  for local acceptor
            labelList donorCellForAcceptor_;

            //- Donor processor for local acceptor
            labelList donorProcForAcceptor_;

            //- Donor processor donor index for local acceptor
            labelList donorProcIndexForAcceptor_;


        // Parallel communication

            //- Map-distribute comms tool
            mutable mapDistribute* mapPtr_;

            //- Local donor cells
            mutable labelList donorCells_;

            //- Donor cell processor (to which donor is sent)
            mutable labelList donorCellsProc_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        oversetAMGInterface(const oversetAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const oversetAMGInterface&);


        //- Init map
        void initMap() const;


    // Private static data

        //- Processor cluster offset index
        static const long procOffset = 12000000;


public:

    //- Runtime type information
    TypeName("overset");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        oversetAMGInterface
        (
            const lduPrimitiveMesh& lduMesh,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    //- Destructor
    virtual ~oversetAMGInterface();


    // Member Functions

        // Access

            //- Return true if interface is coupled
            virtual bool coupled() const
            {
                return true;
            }


        // Communications support


            //- Return the values of the given internal data adjacent to
            //  the interface as a field
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;


        // Agglomeration

            //- Agglomerating given fine-level coefficients
            //  In overset, special treatment is required for upper,
            //  which needs a virtual function.  HJ, 16/Sep/2019
            virtual tmp<scalarField> agglomerateCoeffs
            (
                const scalarField& fineCoeffs
            ) const;


        // Interface transfer functions

            //- Return lduMesh for owner-neighbour addressing
            virtual const lduMesh& ldu() const;

            //- Return access to overset mesh.
            //  Available only on the finest level
            virtual const oversetMesh& overset() const;

            //- Return acceptor cells
            virtual const labelList& acceptorCells() const
            {
                return faceCells_;
            }

            //- Return donor cells for coarse-level addressing
            virtual const labelList& donorCells() const
            {
                return donorCells_;
            }

            //- Return donor cells processor (to which donor is sent)
            virtual const labelList& donorCellsProc() const
            {
                return donorCellsProc_;
            }

            //- Return mapDistribute
            virtual const mapDistribute& map() const;


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

            //- Initialise transfer of internal field adjacent to the interface
            //  For overset, this is fringe cell data
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            //  For overset, this is fringe cell data
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;

            //- Initialise transfer of internal field adjacent to the interface
            //  For overset, this is fringe cell data
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const scalarField& iF
            ) const;

            //- Initialise transfer of internal field adjacent to the interface
            //  For overset, this is fringe cell data
            virtual tmp<scalarField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const scalarField& iF
            ) const;


        //- Overset interface functions

            //- Does this side own the patch ?
            virtual bool master() const;

            //- Is this the fine level?
            virtual bool fineLevel() const;

            //- Return reference tor the oversetLduInterface from which this is
            //  agglomerated
            const oversetLduInterface& fineOversetInterface() const
            {
                return fineOversetInterface_;
            }

            //- Return interface size.  On the coarse level, this is the size
            //  of donor list
            virtual label interfaceSize() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
// #   include "oversetAMGInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
