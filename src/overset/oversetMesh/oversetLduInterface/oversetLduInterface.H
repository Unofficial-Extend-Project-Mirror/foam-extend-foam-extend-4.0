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
    oversetLduInterface

Description
    An abstract base class for overset coupled interfaces.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    oversetLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef oversetLduInterface_H
#define oversetLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class lduMesh;
class oversetMesh;
class mapDistribute;

/*---------------------------------------------------------------------------*\
                      Class oversetLduInterface Declaration
\*---------------------------------------------------------------------------*/

class oversetLduInterface
{

public:

    //- Runtime type information
    TypeName("oversetLduInterface");


    // Constructors

        //- Construct null
        oversetLduInterface()
        {}


    //- Destructor
    virtual ~oversetLduInterface();


    // Member Functions

        //- Is this the master side?
        virtual bool master() const = 0;

        //- Is this the fine level?
        virtual bool fineLevel() const = 0;

        //- Return interface size.  Different on fine and coarse level
        //  due to different comms pattern
        virtual label interfaceSize() const = 0;

        //- Return lduMesh for owner-neighbour addressing
        virtual const lduMesh& ldu() const = 0;

        //- Return oversetMesh for fine-level addressing
        virtual const oversetMesh& overset() const = 0;

        //- Return acceptor cells
        virtual const labelList& acceptorCells() const = 0;

        //- Return donor cells for coarse-level addressing
        virtual const labelList& donorCells() const = 0;

        //- Return donor cells processor (to which donor is sent)
        virtual const labelList& donorCellsProc() const = 0;

        //- Return mapDistribute
        virtual const mapDistribute& map() const = 0;


        // Interface transfer functions

            //- Initialise transfer of internal field adjacent to the interface
            //  For overset, this is fringe cell data
            virtual void initInternalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& iF
            ) const = 0;

            //- Transfer and return internal field adjacent to the interface
            //  For overset, this is fringe cell data
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
