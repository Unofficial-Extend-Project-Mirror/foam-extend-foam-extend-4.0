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
    Foam::cyclicAMGInterfaceField

Description
    AMG agglomerated cyclic interface field.

SourceFiles
    cyclicAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef cyclicAMGInterfaceField_H
#define cyclicAMGInterfaceField_H

#include "AMGInterfaceField.H"
#include "cyclicAMGInterface.H"
#include "cyclicLduInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class cyclicAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class cyclicAMGInterfaceField
:
    public AMGInterfaceField,
    virtual public cyclicLduInterfaceField
{
    // Private data

        //- Local reference cast into the cyclic interface
        const cyclicAMGInterface& cyclicInterface_;

        //- Is the transform required
        bool doTransform_;

        //- Rank of component for transformation
        int rank_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        cyclicAMGInterfaceField(const cyclicAMGInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const cyclicAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("cyclic");


    // Constructors

        //- Construct from AMG interface and fine level interface field
        cyclicAMGInterfaceField
        (
            const AMGInterface& AMGCp,
            const lduInterfaceField& fineInterfaceField
        );


    // Destructor

        virtual ~cyclicAMGInterfaceField();


    // Member Functions

        // Access

            //- Return size
            label size() const
            {
                return cyclicInterface_.size();
            }


        // Coupled interface matrix update

            //- Transform neighbour field
            virtual void transformCoupleField
            (
                scalarField& pnf,
                const direction cmpt
            ) const
            {
                cyclicLduInterfaceField::transformCoupleField(pnf, cmpt);
            }

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


        //- Cyclic interface functions

            //- Does the interface field perform the transfromation
            virtual bool doTransform() const
            {
                return doTransform_;
            }

            //- Return face transformation tensor
            virtual const tensorField& forwardT() const
            {
                return cyclicInterface_.forwardT();
            }

            //- Return neighbour-cell transformation tensor
            virtual const tensorField& reverseT() const
            {
                return cyclicInterface_.reverseT();
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
