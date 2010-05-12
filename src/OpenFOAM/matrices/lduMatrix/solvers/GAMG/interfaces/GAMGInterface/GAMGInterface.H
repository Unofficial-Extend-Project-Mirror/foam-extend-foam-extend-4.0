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
    Foam::GAMGInterface

Description
    Abstract base class for GAMG agglomerated interfaces.

    For complex (interpolative) agglomeration, a single fine-level coefficient
    may be agglomerated into several coarse level coefficients.
    The basic functionality has been rewritten to account for this, using
    three addressing arrays and weights:
    faceCells_: contains coarse interface-to-matrix addressing

    The following three arrays are all of the same size
    fineAddressing_: contains addressing into the fine coefficient array
    restrictAddressing_: contains addressing into the coarse coefficient array
    restrictWeights_: contains restriction factor weights

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
    GAMGInterface.C
    newAmgInterface.C

\*---------------------------------------------------------------------------*/

#ifndef GAMGInterface_H
#define GAMGInterface_H

#include "lduInterface.H"
#include "autoPtr.H"
#include "lduPrimitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class GAMGInterface Declaration
\*---------------------------------------------------------------------------*/

class GAMGInterface
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

        //- Fine addressing
        labelField fineAddressing_;

        //- Restrict addressing
        labelField restrictAddressing_;

        //- Fine level agglomeration weights
        scalarField restrictWeights_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        GAMGInterface(const GAMGInterface&);

        //- Disallow default bitwise assignment
        void operator=(const GAMGInterface&);


public:

    //- Runtime type information
    TypeName("GAMGInterface");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            GAMGInterface,
            lduInterface,
            (
                const lduPrimitiveMesh& lduMesh,
                const lduInterface& fineInterface,
                const labelField& localRestrictAddressing,
                const labelField& neighbourRestrictAddressing
            ),
            (
                lduMesh,
                fineInterface,
                localRestrictAddressing,
                neighbourRestrictAddressing
            )
        );


    // Selectors

        //- Return a pointer to a new interface created on freestore given
        //  the fine interface
        static autoPtr<GAMGInterface> New
        (
            const lduPrimitiveMesh& lduMesh,
            const lduInterface& fineInterface,
            const labelField& localRestrictAddressing,
            const labelField& neighbourRestrictAddressing
        );


    // Constructors

        //- Construct from fine-level interface,
        //  local and neighbour restrict addressing
        GAMGInterface(const lduPrimitiveMesh& lduMesh)
        :
            lduMesh_(lduMesh)
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

            //- Merge the next level with this level
            //  combining the face-restrict addressing
            //  and copying the face-cell addressing
            void combine(const GAMGInterface&);

            //- Agglomerating the given fine-level coefficients and return
            virtual tmp<scalarField> agglomerateCoeffs
            (
                const scalarField& fineCoeffs
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "GAMGInterfaceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
