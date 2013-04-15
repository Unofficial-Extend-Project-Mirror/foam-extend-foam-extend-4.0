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

#ifndef tensorProcessorBlockGAMGInterfaceField_H
#define tensorProcessorBlockGAMGInterfaceField_H

#include "processorBlockGAMGInterfaceField.H"
#include "BlockGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "tensorCoeffField.H"
#include "BlockCoeff.H"
#include "tensorBlockCoeff.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

template<>                                                                     
void processorBlockGAMGInterfaceField<tensor>::initInterfaceMatrixUpdate                     
(                                                                               
    const Field<tensor>& psiInternal,                                             
    Field<tensor>&,                                                               
    const BlockLduMatrix<tensor>&,                                                
    const CoeffField<tensor>&,                                                    
    const Pstream::commsTypes commsType                                         
) const                                                                         
{                                                                               
    FatalErrorIn("Foam::tensor Foam::processorBlockGAMGInterfaceField<tensor>::initInterafaceMatrix(...)")
        << "Not implemented" << abort(FatalError);
}                                                                               
                                                                                
template<>                                                                     
void processorBlockGAMGInterfaceField<tensor>::updateInterfaceMatrix                         
(                                                                               
    const Field<tensor>& psiInternal,                                             
    Field<tensor>& result,                                                        
    const BlockLduMatrix<tensor>&,                                                
    const CoeffField<tensor>& coeffs,                                             
    const Pstream::commsTypes commsType                                         
) const                                                                         
{                                                                               
    FatalErrorIn("Foam::tensor Foam::processorBlockGAMGInterfaceField<tensor>::initInterafaceMatrix(...)")
        << "Not implemented" << abort(FatalError);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
