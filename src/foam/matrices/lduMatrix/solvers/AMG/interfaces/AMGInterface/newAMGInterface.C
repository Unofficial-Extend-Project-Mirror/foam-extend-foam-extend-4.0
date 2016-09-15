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

\*---------------------------------------------------------------------------*/

#include "AMGInterface.H"
#include "lduMatrix.H"


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::AMGInterface> Foam::AMGInterface::New
(
    const lduPrimitiveMesh& lduMesh,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
{
    word coupleType(fineInterface.type());

    lduInterfaceConstructorTable::iterator cstrIter =
        lduInterfaceConstructorTablePtr_->find(coupleType);

    if (cstrIter == lduInterfaceConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "AMGInterface::New\n"
            "(\n"
            "    const lduPrimitiveMesh& lduMesh,\n"
            "    const lduInterfacePtrsList& coarseInterfaces,\n"
            "    const lduInterface& fineInterface,\n"
            "    const labelField& localRestrictAddressing,\n"
            "    const labelField& neighbourRestrictAddressing\n"
            ")"
        )   << "Unknown AMGInterface type " << coupleType << ".\n"
            << "Valid AMGInterface types are :"
            << lduInterfaceConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<AMGInterface>
    (
        cstrIter()
        (
            lduMesh,
            coarseInterfaces,
            fineInterface,
            localRestrictAddressing,
            neighbourRestrictAddressing
        )
    );
}


// ************************************************************************* //
