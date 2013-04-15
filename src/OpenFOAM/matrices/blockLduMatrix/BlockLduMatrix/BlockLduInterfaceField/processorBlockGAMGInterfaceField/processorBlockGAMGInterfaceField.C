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

\*---------------------------------------------------------------------------*/

#include "processorBlockGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::processorBlockGAMGInterfaceField<Type>::processorBlockGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const BlockLduInterfaceField<Type>& fineInterfaceField
)
:
    BlockGAMGInterfaceField<Type>(GAMGCp, fineInterfaceField),
    procInterface_(refCast<const processorGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0)
{
    // If the interface based on a patch this must be taken care specially of
    if (isA<processorBlockLduInterfaceField<Type> >(fineInterfaceField))
    { 
        const processorBlockLduInterfaceField<Type>& p =
            refCast<const processorBlockLduInterfaceField<Type> >(fineInterfaceField);

        doTransform_ = p.doTransform();
        rank_ = p.rank();
    }
    else
    {
        FatalErrorIn("Foam::processorBlockGAMGInterfaceField<Type> Constructor")
            << "fineInterface must be of processor type and either" << endl 
            << "    processorBlockLduInterfaceField<Type> or " << endl
            << "    processorFvPatchField<Type> " << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::processorBlockGAMGInterfaceField<Type>::~processorBlockGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>                                                                     
void Foam::processorBlockGAMGInterfaceField<Type>::initInterfaceMatrixUpdate                     
(                                                                               
    const Field<Type>& psiInternal,                                             
    Field<Type>&,                                                               
    const BlockLduMatrix<Type>&,                                                
    const CoeffField<Type>&,                                                    
    const Pstream::commsTypes commsType                                         
) const                                                                         
{                                                                               
    procInterface_.compressedSend                                                   
    (                                                                           
        commsType,                                                              
        procInterface_.interfaceInternalField(psiInternal)()
    );                                                                          
}                                                                               
                                                                                
template<class Type>                                                                     
void Foam::processorBlockGAMGInterfaceField<Type>::updateInterfaceMatrix                         
(                                                                               
    const Field<Type>& psiInternal,                                             
    Field<Type>& result,                                                        
    const BlockLduMatrix<Type>&,                                                
    const CoeffField<Type>& coeffs,                                             
    const Pstream::commsTypes commsType                                         
) const                                                                         
{                                                                               
    Field<Type> pnf
    (
        coeffs.size()
    );                                                                           
 
    if (coeffs.activeType() == blockCoeffBase::SCALAR)                          
    {                                                                           
        pnf = coeffs.asScalar() *                                               
            procInterface_.compressedReceive<Type>(commsType, this->size())();      
    }                                                                           
    else if (coeffs.activeType() == blockCoeffBase::LINEAR)                     
    {                                                                           
        pnf = cmptMultiply(coeffs.asLinear(),                                   
            procInterface_.compressedReceive<Type>(commsType, this->size())()       
        );                                                                      
    }                                                                           
    else if (coeffs.activeType() == blockCoeffBase::SQUARE)                     
    {                                                                           
        pnf = coeffs.asSquare() &                                               
            procInterface_.compressedReceive<Type>(commsType, this->size())();      
    }                                                                           
                                                                                
    const unallocLabelList& faceCells = procInterface_.faceCells();              
    
    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= pnf[elemI];
    }                                                                                
}


// ************************************************************************* //
