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
    Foam::lduInterfaceField

Description
    An abstract base class for implicitly-coupled interface fields
    e.g. processor and cyclic patch fields.

SourceFiles
    lduInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef lduInterfaceField_H
#define lduInterfaceField_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class lduMatrix;

/*---------------------------------------------------------------------------*\
                     Class lduInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class lduInterfaceField
{
    // Private data

        //- Reference to the coupled patch this field is defined for
        const lduInterface& coupledInterface_;

        //- Update index used so that updateInterfaceMatrix is called only once
        //  during the construction of the matrix
        bool updatedMatrix_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        lduInterfaceField(const lduInterfaceField&);

        //- Disallow default bitwise assignment
        void operator=(const lduInterfaceField&);


public:

    //- Runtime type information
    TypeName("lduInterfaceField");


    // Constructors

        //- Construct given coupled patch
        lduInterfaceField(const lduInterface& patch)
        :
            coupledInterface_(patch),
            updatedMatrix_(false)
        {}


    //- Destructor
    virtual ~lduInterfaceField();


    // Member Functions

        // Access

            //- Return the interface
            const lduInterface& coupledInterface() const
            {
                return coupledInterface_;
            }

            //- Return the interface type
            virtual const word& interfaceFieldType() const
            {
                return type();
            }


        // Coupled interface matrix update

            //- Whether matrix has been updated
            bool updatedMatrix() const
            {
                return updatedMatrix_;
            }

            //- Whether matrix has been updated
            bool& updatedMatrix()
            {
                return updatedMatrix_;
            }

            //- Is all data available
            virtual bool ready() const
            {
                return true;
            }

            //- Transform neighbour field
            virtual void transformCoupleField
            (
                scalarField& pnf,
                const direction cmpt
            ) const = 0;

            //- Initialise neighbour matrix update
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const
            {}

            //- Update result field based on interface functionality
            virtual void updateInterfaceMatrix
            (
                const scalarField&,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const Pstream::commsTypes commsType,
                const bool switchToLhs
            ) const = 0;

            //- Update the constraints
            virtual void updateConstraints
            (
                const scalarField&,
                scalarField&,
                const label
            ) const
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
