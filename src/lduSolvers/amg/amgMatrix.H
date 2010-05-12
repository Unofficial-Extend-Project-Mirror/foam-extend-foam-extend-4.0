/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    amgMatrix

Description
    Matrix wrapper class containing ldu addressing, matrix, coupling interfaces
    coupling interface fields and coupling coefficients

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    amgMatrixI.H
    amgMatrix.C
    amgMatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef amgMatrix_H
#define amgMatrix_H

#include "lduPrimitiveMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class amgMatrix Declaration
\*---------------------------------------------------------------------------*/

class amgMatrix
{
    // Private data

        //- Matrix addressing object
        lduPrimitiveMesh* addrPtr_;

        //- Coupling interfaces fields
        lduInterfacePtrsList* interfacesPtr_;

         //- LDU matrix
        lduMatrix* matrixPtr_;

        //- Coupling interface fields
        lduInterfaceFieldPtrsList* interfaceFieldsPtr_;

        //- Coupling coefficients, upper
        FieldField<Field, scalar>* coupleBouCoeffsPtr_;

        //- Coupling coefficients, lower
        FieldField<Field, scalar>* coupleIntCoeffsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        amgMatrix(const amgMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const amgMatrix&);


public:

    // Static data members


    // Constructors

        //- Construct from components
        amgMatrix
        (
            lduPrimitiveMesh* addrPtr,
            lduInterfacePtrsList* interfacesPtr,
            lduMatrix* matrixPtr,
            FieldField<Field, scalar>* coupleBouCoeffsPtr,
            FieldField<Field, scalar>* coupleIntCoeffsPtr,
            lduInterfaceFieldPtrsList* interfaceFieldsPtr
        )
        :
            addrPtr_(addrPtr),
            interfacesPtr_(interfacesPtr),
            matrixPtr_(matrixPtr),
            interfaceFieldsPtr_(interfaceFieldsPtr),
            coupleBouCoeffsPtr_(coupleBouCoeffsPtr),
            coupleIntCoeffsPtr_(coupleIntCoeffsPtr)
        {}


    // Destructor

        ~amgMatrix()
        {
            deleteDemandDrivenData(addrPtr_);

            // Clear the interface storage by hand.  It is a list of ptrs
            // not a PtrList for consistency of the interface
            if (interfacesPtr_)
            {
                lduInterfacePtrsList& curLevel = *interfacesPtr_;

                forAll (curLevel, i)
                {
                    if (curLevel.set(i))
                    {
                        delete curLevel(i);
                    }
                }
            }

            deleteDemandDrivenData(interfacesPtr_);

            deleteDemandDrivenData(matrixPtr_);

            // Clear the interface field storage by hand.  It is a list of ptrs
            // not a PtrList for consistency of the interface
            if (interfaceFieldsPtr_)
            {
                lduInterfaceFieldPtrsList& curLevelField = *interfaceFieldsPtr_;

                forAll (curLevelField, i)
                {
                    if (curLevelField.set(i))
                    {
                        delete curLevelField(i);
                    }
                }
            }

            deleteDemandDrivenData(interfaceFieldsPtr_);
            deleteDemandDrivenData(coupleBouCoeffsPtr_);
            deleteDemandDrivenData(coupleIntCoeffsPtr_);
        }


    // Member Functions

        //- Return size
        label size() const
        {
            return addrPtr_->size();
        }

        //- Return matrix
        const lduMatrix& matrix() const
        {
            return *matrixPtr_;
        }

        //- Return coupling interfaces
        const lduInterfacePtrsList& interfaces() const
        {
            return *interfacesPtr_;
        }

        //- Return coupling interfaceFields
        const lduInterfaceFieldPtrsList& interfaceFields() const
        {
            return *interfaceFieldsPtr_;
        }

        //- Return boundary coefficients
        const FieldField<Field, scalar>& coupleBouCoeffs() const
        {
            return *coupleBouCoeffsPtr_;
        }

        //- Return internal coefficients
        const FieldField<Field, scalar>& coupleIntCoeffs() const
        {
            return *coupleIntCoeffsPtr_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
