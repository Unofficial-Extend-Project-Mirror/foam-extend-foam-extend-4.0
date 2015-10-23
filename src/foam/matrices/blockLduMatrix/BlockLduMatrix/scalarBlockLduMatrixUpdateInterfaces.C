/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description
    Update of block interfaces

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TODO - This code is currently not called so we have specialized
// initInterfaceMatrixUpdate in processorFvPatchScalarfield
// This needs to be fixed

template<>
void Foam::BlockLduMatrix<scalar>::initInterfaces
(
    const FieldField<CoeffField, scalar>& coupleCoeffs,
    scalarField& result,
    const scalarField& psi,
    const bool switchToLhs
) const
{
    if
    (
        Pstream::defaultCommsType() == Pstream::blocking
     || Pstream::defaultCommsType() == Pstream::nonBlocking
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
                    Pstream::defaultCommsType,
                    switchToLhs
                );
            }
        }
    }
    else if (Pstream::defaultCommsType() == Pstream::scheduled)
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
                    Pstream::blocking,
                    switchToLhs
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
    const scalarField& psi,
    const bool switchToLhs
) const
{
    if
    (
        Pstream::defaultCommsType() == Pstream::blocking
     || Pstream::defaultCommsType() == Pstream::nonBlocking
    )
    {
        // Block until all sends/receives have been finished
        if (Pstream::defaultCommsType() == Pstream::nonBlocking)
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
                    Pstream::defaultCommsType,
                    switchToLhs
                );
            }
        }
    }
    else if (Pstream::defaultCommsType() == Pstream::scheduled)
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
                        Pstream::scheduled,
                        switchToLhs
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
                        Pstream::scheduled,
                        switchToLhs
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
                    Pstream::blocking,
                    switchToLhs
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
