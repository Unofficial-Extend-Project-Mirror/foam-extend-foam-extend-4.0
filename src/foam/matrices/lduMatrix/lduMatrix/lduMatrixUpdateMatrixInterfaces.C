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

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduMatrix::initMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt,
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
                    psiif,
                    result,
                    *this,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    static_cast<Pstream::commsTypes>
                    (
                        Pstream::defaultCommsType()
                    ),
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
            label interfaceI = patchSchedule.size()/2;
            interfaceI < interfaces.size();
            interfaceI++
        )
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].initInterfaceMatrixUpdate
                (
                    psiif,
                    result,
                    *this,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    Pstream::blocking,
                    switchToLhs
                );
            }
        }
    }
    else
    {
        FatalErrorIn("lduMatrix::initMatrixInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType()]
            << exit(FatalError);
    }
}


void Foam::lduMatrix::updateMatrixInterfaces
(
    const FieldField<Field, scalar>& coupleCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const scalarField& psiif,
    scalarField& result,
    const direction cmpt,
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
                    psiif,
                    result,
                    *this,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    static_cast<Pstream::commsTypes>
                    (
                        Pstream::defaultCommsType()
                    ),
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
                        psiif,
                        result,
                        *this,
                        coupleCoeffs[interfaceI],
                        cmpt,
                        Pstream::scheduled,
                        switchToLhs
                    );
                }
                else
                {
                    interfaces[interfaceI].updateInterfaceMatrix
                    (
                        psiif,
                        result,
                        *this,
                        coupleCoeffs[interfaceI],
                        cmpt,
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
            interfaceI < interfaces.size();
            interfaceI++
        )
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].updateInterfaceMatrix
                (
                    psiif,
                    result,
                    *this,
                    coupleCoeffs[interfaceI],
                    cmpt,
                    Pstream::blocking,
                    switchToLhs
                );
            }
        }
    }
    else
    {
        FatalErrorIn("lduMatrix::updateMatrixInterfaces")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType()]
            << exit(FatalError);
    }
}


// ************************************************************************* //
