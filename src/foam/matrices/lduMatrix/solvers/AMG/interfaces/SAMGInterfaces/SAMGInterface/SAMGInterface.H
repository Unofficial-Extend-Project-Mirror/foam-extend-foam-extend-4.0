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
    Foam::SAMGInterface

Description
    Abstract base class for AMG selected interfaces.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    SAMGInterface.C
    newAmgInterface.C

\*---------------------------------------------------------------------------*/

#ifndef SAMGInterface_H
#define SAMGInterface_H

#include "lduInterface.H"
#include "autoPtr.H"
#include "lduPrimitiveMesh.H"
#include "coeffFields.H"

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


protected:

    // Protected data

        //- Face-cell addressing.  Contains coarse level addressing
        labelField faceCells_;

        //- Fine addressing.  Contains fine index for each coarse face
        labelField fineAddressing_;


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
                const lduInterfacePtrsList& coarseInterfaces,
                const lduInterface& fineInterface,
                const labelField& localSelectAddressing,
                const labelField& neighbourSelectAddressing
            ),
            (
                lduMesh,
                coarseInterfaces,
                fineInterface,
                localSelectAddressing,
                neighbourSelectAddressing
            )
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<SAMGInterface> New
        (
            const lduPrimitiveMesh& lduMesh,
            const lduInterfacePtrsList& coarseInterfaces,
            const lduInterface& fineInterface,
            const labelField& localSelectAddressing,
            const labelField& neighbourSelectAddressing
        );


    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour select addressing
        SAMGInterface(const lduPrimitiveMesh& lduMesh)
        :
            lduMesh_(lduMesh)
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

            //- Return size
            virtual label size() const
            {
                return faceCells_.size();
            }

            //- Return faceCell addressing
            virtual const unallocLabelList& faceCells() const
            {
                return faceCells_;
            }

            //- Return fine addressing
            const labelField& fineAddressing() const
            {
                return fineAddressing_;
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


        // Agglomeration

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
