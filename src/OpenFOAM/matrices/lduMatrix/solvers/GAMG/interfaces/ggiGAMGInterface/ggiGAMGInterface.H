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
    Foam::ggiGAMGInterface

Description
    GAMG agglomerated ggi interface.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    ggiGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef ggiGAMGInterface_H
#define ggiGAMGInterface_H

#include "GAMGInterface.H"
#include "ggiLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class ggiGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class ggiGAMGInterface
:
    public GAMGInterface,
    virtual public ggiLduInterface
{
    // Private data

        //- Reference tor the ggiLduInterface from which this is
        //  agglomerated
        const ggiLduInterface& fineGgiInterface_;

        //- Zone size
        label zoneSize_;

        //- Zone addressing
        labelList zoneAddressing_;


        // Parallel communication optimisation, stored on master processor

            //- Has reduce been initialised?
            mutable bool initReduce_;

            //- List of zone faces indices received from each processor
            mutable labelListList receiveAddr_;

            //- List of zone faces indices to send to each processor
            mutable labelListList sendAddr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ggiGAMGInterface(const ggiGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const ggiGAMGInterface&);


        //- Init fast reduce
        void initFastReduce() const;

        //- Fast reduce
        //  Note: contains global communications
        //  New, HJ, 24/Jun/2011
        template<class Type>
        tmp<Field<Type> > fastReduce(const UList<Type>&) const;

        //- Expand data to zone size
        //  Note: contains global reduce.  HJ, 15/May/2009
        //  Obsolete.  HJ, 24/Jun/2011
        template<class Type>
        tmp<Field<Type> > expand(const UList<Type>&) const;


    // Private static data

        //- Processor cluster offset index
        static const label procOffset = 1000000;


public:

    //- Runtime type information
    TypeName("ggi");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        ggiGAMGInterface
        (
            const lduPrimitiveMesh& lduMesh,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    // Destructor

        virtual ~ggiGAMGInterface();


    // Member Functions

        // Agglomeration

            //- Agglomerating the given fine-level coefficients and return
            virtual tmp<scalarField> agglomerateCoeffs
            (
                const scalarField& fineCoeffs
            ) const;


        // Interface transfer functions

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

            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const scalarField& iF
            ) const;

            //- Transfer and return internal field adjacent to the interface
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const;

            virtual tmp<scalarField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const scalarField& iF
            ) const;


        //- GGI interface functions

            //- Is this the master side?
            virtual bool master() const;

            //- Is this the fine level?
            virtual bool fineLevel() const;

            //- Return shadow patch index
            virtual label shadowIndex() const;

            //- Return shadow interface
            virtual const ggiLduInterface& shadowInterface() const;

            //- Return zone size
            virtual label zoneSize() const;

            //- Return zone addressing
            virtual const labelList& zoneAddressing() const;

            //- Return addressing
            virtual const labelListList& addressing() const;

            //- Is the patch localised on a single processor
            virtual bool localParallel() const;

            //- Return weights
            virtual const scalarListList& weights() const;

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const;

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ggiGAMGInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
