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
    Foam::mixingPlaneGAMGInterface

Description
    GAMG agglomerated mixingPlane interface.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlaneGAMGInterface.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneGAMGInterface_H
#define mixingPlaneGAMGInterface_H

#include "GAMGInterface.H"
#include "mixingPlaneLduInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class mixingPlaneGAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class mixingPlaneGAMGInterface
:
    public GAMGInterface,
    virtual public mixingPlaneLduInterface
{
    // Private data

        //- Reference tor the mixingPlaneLduInterface from which this is
        //  agglomerated
        const mixingPlaneLduInterface& fineMixingPlaneInterface_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixingPlaneGAMGInterface(const mixingPlaneGAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const mixingPlaneGAMGInterface&);


    // Private static data

        //- Processor cluster offset index
        static const label procOffset = 1000000;


public:

    //- Runtime type information
    TypeName("mixingPlane");


    // Constructors

        //- Construct from fine level interface,
        //  local and neighbour restrict addressing
        mixingPlaneGAMGInterface
        (
            const lduPrimitiveMesh& lduMesh,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    // Destructor

        virtual ~mixingPlaneGAMGInterface();


    // Member Functions

        // Access

            //- Return true if interface is coupled
            virtual bool coupled() const
            {
                return true;
            }


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


        //- MixingPlane interface functions

            //- Is this the master side?
            virtual bool master() const;

            //- Return shadow patch index
            virtual label shadowIndex() const;

            //- Return shadow interface
            virtual const mixingPlaneLduInterface& shadowInterface() const;

            //- Return addressing
            virtual const labelListList& addressing() const;

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

#endif

// ************************************************************************* //
