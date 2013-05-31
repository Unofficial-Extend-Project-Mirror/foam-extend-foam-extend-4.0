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
    Foam::mixingPlaneGAMGInterfaceField

Description
    GAMG agglomerated mixingPlane interface field.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlaneGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneGAMGInterfaceField_H
#define mixingPlaneGAMGInterfaceField_H

#include "GAMGInterfaceField.H"
#include "mixingPlaneGAMGInterface.H"
#include "mixingPlaneLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class mixingPlaneGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class mixingPlaneGAMGInterfaceField
:
    public GAMGInterfaceField,
    virtual public mixingPlaneLduInterfaceField
{
    // Private data

        //- Local reference cast into the mixingPlane interface
        const mixingPlaneGAMGInterface& mixingPlaneInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;

        //- Transfer buffer
        mutable scalarField transferBuffer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        mixingPlaneGAMGInterfaceField(const mixingPlaneGAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const mixingPlaneGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("mixingPlane");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        mixingPlaneGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterfaceField
        );


    // Destructor

        virtual ~mixingPlaneGAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return mixingPlaneInterface_.size();
            }


        // Interface matrix update

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix& m,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const;


        //- MixingPlane interface functions

            //- Does the interface field perform the transfromation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return mixingPlaneInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return mixingPlaneInterface_.reverseT();
            }

            //- Return rank of component for transform
            virtual int rank() const
            {
                return rank_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
