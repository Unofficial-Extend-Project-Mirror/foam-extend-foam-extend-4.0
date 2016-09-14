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
    mixingPlaneLduInterface

Description
    An abstract base class for mixingPlane coupled interfaces.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlaneLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneLduInterface_H
#define mixingPlaneLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class mixingPlaneLduInterface Declaration
\*---------------------------------------------------------------------------*/

class mixingPlaneLduInterface
{
protected:

    // Protected data

        //- Transfer buffer
        mutable labelField labelTransferBuffer_;

        //- Field transfer buffer
        mutable scalarField fieldTransferBuffer_;


public:

    //- Runtime type information
    TypeName("mixingPlaneLduInterface");


    // Constructors

        //- Construct null
        mixingPlaneLduInterface()
        :
            labelTransferBuffer_(),
            fieldTransferBuffer_()
        {}


    //- Destructor
    virtual ~mixingPlaneLduInterface();


    // Member Functions

        //- Is this the master side?
        virtual bool master() const = 0;

        //- Return shadow interface index
        virtual label shadowIndex() const = 0;

        //- Return shadow interface
        virtual const mixingPlaneLduInterface& shadowInterface() const = 0;

        //- Return fine level addressing
        virtual const labelListList& addressing() const = 0;

        //- Return weights
        virtual const scalarListList& weights() const = 0;

        //- Return face transformation tensor
        virtual const tensorField& forwardT() const = 0;

        //- Return face reverse transformation tensor
        virtual const tensorField& reverseT() const = 0;


        // Communications support

            //- Return communicator used for parallel communication
            virtual int comm() const = 0;

            //- Return message tag used for sending
            virtual int tag() const = 0;


        // Transfer buffer access

            //- Return contents of the label transfer buffer
            const labelField& labelTransferBuffer() const
            {
                return labelTransferBuffer_;
            }

            //- Return contents of the field transfer buffer
            const scalarField& fieldTransferBuffer() const
            {
                return fieldTransferBuffer_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
