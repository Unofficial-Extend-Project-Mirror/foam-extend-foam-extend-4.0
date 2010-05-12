// The FOAM Project // File: updateCoupledInterfaces.C
/*
-------------------------------------------------------------------------------
 =========         | Class Implementation
 \\      /         |
  \\    /          | Name:   updateCoupledInterfaces
   \\  /           | Family: matrix
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION

AUTHOR
    Hrvoje Jasak and Henry G. Weller.

-------------------------------------------------------------------------------
*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void lduMatrix::updateMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduCoupledInterfacePtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt
) const
{
    // Initialise coupling
    forAll (interfaces, interfaceI)
    {
        if (interfaces[interfaceI]->coupled())
        {
            interfaces[interfaceI]->initInterfaceMatrixUpdate
            (
                psiif,
                result,
                *this,
                coupleCoeffs[interfaceI],
                cmpt,
                false
            );
        }
    }

    // Update coupled interface
    forAll (interfaces, interfaceI)
    {
        if (interfaces[interfaceI]->coupled())
        {
            interfaces[interfaceI]->updateInterfaceMatrix
            (
                psiif,
                result,
                *this,
                coupleCoeffs[interfaceI],
                cmpt
            );
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
