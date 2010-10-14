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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Update of block interfaces

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO - This code is currently not called so we have specialized
// initInterfaceMatrixUpdate in processorFvPatchScalarfield. This needs to be fixed

template<>
void Foam::BlockLduMatrix<scalar>::initInterfaces
(
    const FieldField<CoeffField, scalar>& coupleCoeffs,
    scalarField& result,
    const scalarField& psi
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        forAll (interfaces, interfaceI)
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].initInterfaceMatrixUpdate
                (
                    psi,
                    result,
                    *this,
                    coupleCoeffs[interfaceI].asScalar(),
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfaceI=patchSchedule.size()/2;
            interfaceI<interfaces.size();
            interfaceI++
        )
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].initInterfaceMatrixUpdate
                (
                    psi,
                    result,
                    *this,
                    coupleCoeffs[interfaceI].asScalar(),
                    Pstream::blocking
                );
            }
        }
    }    
    else
    {
        FatalErrorIn("BlockLduMatrix<scalar>::initMatrixInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


template<>
void Foam::BlockLduMatrix<scalar>::updateInterfaces
(
    const FieldField<CoeffField, scalar>& coupleCoeffs,
    scalarField& result,
    const scalarField& psi
) const
{
    if
    (
        Pstream::defaultCommsType == Pstream::blocking
     || Pstream::defaultCommsType == Pstream::nonBlocking
    )
    {
        // Block until all sends/receives have been finished
        if (Pstream::defaultCommsType == Pstream::nonBlocking)
        {
            IPstream::waitRequests();
            OPstream::waitRequests();
        }

        forAll (interfaces, interfaceI)
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].updateInterfaceMatrix
                (
                    psi,
                    result,
                    *this,
                    coupleCoeffs[interfaceI].asScalar(),
                    Pstream::defaultCommsType
                );
            }
        }
    }
    else if (Pstream::defaultCommsType == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        forAll (patchSchedule, i)
        {
            label interfaceI = patchSchedule[i].patch;

            if (interfaces.set(interfaceI))
            {
                if (patchSchedule[i].init)
                {
                    interfaces[interfaceI].initInterfaceMatrixUpdate
                    (
                        psi,
                        result,
                        *this,
                        coupleCoeffs[interfaceI].asScalar(),
                        Pstream::scheduled
                    );
                }
                else
                {
                    interfaces[interfaceI].updateInterfaceMatrix
                    (
                        psi,
                        result,
                        *this,
                        coupleCoeffs[interfaceI].asScalar(),
                        Pstream::scheduled
                    );
                }
            }
        }

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfaceI=patchSchedule.size()/2;
            interfaceI<interfaces.size();
            interfaceI++
        )
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].updateInterfaceMatrix
                (
                    psi,
                    result,
                    *this,
                    coupleCoeffs[interfaceI].asScalar(),
                    Pstream::blocking
                );
            }
        }
    }
    else
    {
        FatalErrorIn("BlockLduMatrix<scalar>::updateMatrixInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


// ************************************************************************* //
