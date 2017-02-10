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

Description
    Update of block interfaces

\*---------------------------------------------------------------------------*/

#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::BlockLduMatrix<Type>::initInterfaces
(
    const FieldField<CoeffField, Type>& interfaceCoeffs,
    TypeField& result,
    const TypeField& psi,
    const bool switchToLhs
) const
{
    if
    (
        Pstream::defaultComms() == Pstream::blocking
     || Pstream::defaultComms() == Pstream::nonBlocking
    )
    {
        forAll (interfaces_, interfaceI)
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].initInterfaceMatrixUpdate
                (
                    psi,
                    result,
                    *this,
                    interfaceCoeffs[interfaceI],
                    Pstream::defaultComms(),
                    switchToLhs
                );
            }
        }
    }
    else if (Pstream::defaultComms() == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over the "global" patches are on the list of interfaces but
        // beyond the end of the schedule which only handles "normal" patches
        for
        (
            label interfaceI = patchSchedule.size()/2;
            interfaceI < interfaces_.size();
            interfaceI++
        )
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].initInterfaceMatrixUpdate
                (
                    psi,
                    result,
                    *this,
                    interfaceCoeffs[interfaceI],
                    Pstream::blocking,
                    switchToLhs
                );
            }
        }
    }
    else
    {
        FatalErrorIn("BlockLduMatrix<Type>::initMatrixInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType()]
            << exit(FatalError);
    }
}


template<class Type>
void Foam::BlockLduMatrix<Type>::updateInterfaces
(
    const FieldField<CoeffField, Type>& interfaceCoeffs,
    TypeField& result,
    const TypeField& psi,
    const bool switchToLhs
) const
{
    if (Pstream::defaultComms() == Pstream::blocking)
    {
        forAll (interfaces_, interfaceI)
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].updateInterfaceMatrix
                (
                    psi,
                    result,
                    *this,
                    interfaceCoeffs[interfaceI],
                    Pstream::defaultComms(),
                    switchToLhs
                );
            }
        }
    }
    else if (Pstream::defaultComms() == Pstream::nonBlocking)
    {
        // Try and consume interfaces as they become available
        bool allUpdated = false;

        for (label i = 0; i < Pstream::nPollProcInterfaces; i++)
        {
            allUpdated = true;

            forAll (interfaces_, interfaceI)
            {
                if (interfaces_.set(interfaceI))
                {
                    if (!interfaces_[interfaceI].updatedMatrix())
                    {
                        if (interfaces_[interfaceI].ready())
                        {
                            interfaces_[interfaceI].updateInterfaceMatrix
                            (
                                psi,
                                result,
                                *this,
                                interfaceCoeffs[interfaceI],
                                Pstream::defaultComms(),
                                switchToLhs
                            );
                        }
                        else
                        {
                            allUpdated = false;
                        }
                    }
                }
            }

            if (allUpdated)
            {
                break;
            }
        }

        // Block for everything
        if (Pstream::parRun())
        {
            if (allUpdated)
            {
                // All received. Just remove all storage of requests
                // Note that we don't know what starting number of requests
                // was before start of sends and receives (since set from
                // initMatrixInterfaces) so set to 0 and lose any in-flight
                // requests.
                Pstream::resetRequests(0);
            }
            else
            {
                // Block for all requests and remove storage
                IPstream::waitRequests();
                OPstream::waitRequests();
            }
        }

        // Consume
        forAll (interfaces_, interfaceI)
        {
            if
            (
                interfaces_.set(interfaceI)
            && !interfaces_[interfaceI].updatedMatrix()
            )
            {
                interfaces_[interfaceI].updateInterfaceMatrix
                (
                    psi,
                    result,
                    *this,
                    interfaceCoeffs[interfaceI],
                    Pstream::defaultComms(),
                    switchToLhs
                );
            }
        }
    }
    else if (Pstream::defaultComms() == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule = this->patchSchedule();

        // Loop over all the "normal" interfaces relating to standard patches
        forAll (patchSchedule, i)
        {
            label interfaceI = patchSchedule[i].patch;

            if (interfaces_.set(interfaceI))
            {
                if (patchSchedule[i].init)
                {
                    interfaces_[interfaceI].initInterfaceMatrixUpdate
                    (
                        psi,
                        result,
                        *this,
                        interfaceCoeffs[interfaceI],
                        Pstream::scheduled,
                        switchToLhs
                    );
                }
                else
                {
                    interfaces_[interfaceI].updateInterfaceMatrix
                    (
                        psi,
                        result,
                        *this,
                        interfaceCoeffs[interfaceI],
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
            label interfaceI = patchSchedule.size()/2;
            interfaceI < interfaces_.size();
            interfaceI++
        )
        {
            if (interfaces_.set(interfaceI))
            {
                interfaces_[interfaceI].updateInterfaceMatrix
                (
                    psi,
                    result,
                    *this,
                    interfaceCoeffs[interfaceI],
                    Pstream::blocking,
                    switchToLhs
                );
            }
        }
    }
    else
    {
        FatalErrorIn("BlockLduMatrix<Type>::updateInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType()]
            << exit(FatalError);
    }
}


// ************************************************************************* //
