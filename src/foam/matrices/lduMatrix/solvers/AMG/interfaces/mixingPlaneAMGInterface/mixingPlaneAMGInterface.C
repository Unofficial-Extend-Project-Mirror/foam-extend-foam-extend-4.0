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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, 2009.

\*---------------------------------------------------------------------------*/

#include "mixingPlaneAMGInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingPlaneAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        AMGInterface,
        mixingPlaneAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingPlaneAMGInterface::mixingPlaneAMGInterface
(
    const lduPrimitiveMesh& lduMesh,
    const lduInterfacePtrsList& coarseInterfaces,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
:
    AMGInterface(lduMesh),
    fineMixingPlaneInterface_
    (
        refCast<const mixingPlaneLduInterface>(fineInterface)
    ),
    comm_(fineMixingPlaneInterface_.comm()),
    tag_(fineMixingPlaneInterface_.tag())
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::mixingPlaneAMGInterface::~mixingPlaneAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::mixingPlaneAMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    // AMG agglomeration missing
    notImplemented("mixingPlaneAMGInterface::agglomerateCoeffs");
    tmp<scalarField> tcoarseCoeffs(new scalarField(size(), 0.0));

    return tcoarseCoeffs;
}


bool Foam::mixingPlaneAMGInterface::master() const
{
    return fineMixingPlaneInterface_.master();
}


Foam::label Foam::mixingPlaneAMGInterface::shadowIndex() const
{
    return fineMixingPlaneInterface_.shadowIndex();
}


const Foam::mixingPlaneLduInterface&
Foam::mixingPlaneAMGInterface::shadowInterface() const
{
    return refCast<const mixingPlaneLduInterface>
    (
        ldu().interfaces()[shadowIndex()]
    );
}


const Foam::labelListList& Foam::mixingPlaneAMGInterface::addressing() const
{
    FatalErrorIn
    (
        "const labelListList& mixingPlaneAMGInterface::addressing() const"
    )   << "Requested fine addressing at coarse level"
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::mixingPlaneAMGInterface::weights() const
{
    FatalErrorIn
    (
        "const labelListList& mixingPlaneAMGInterface::weights() const"
    )   << "Requested fine addressing at coarse level"
        << abort(FatalError);
    return scalarListList::null();
}


const Foam::tensorField& Foam::mixingPlaneAMGInterface::forwardT() const
{
    return fineMixingPlaneInterface_.forwardT();
}


const Foam::tensorField& Foam::mixingPlaneAMGInterface::reverseT() const
{
    return fineMixingPlaneInterface_.reverseT();
}


void Foam::mixingPlaneAMGInterface::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::mixingPlaneAMGInterface::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    return this->shadowInterface().labelTransferBuffer();
}


void Foam::mixingPlaneAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    // NOTE: Change this: requires fast reduce.  HJ, 13/Jun/20106
    labelTransferBuffer_ = interfaceInternalField(iF);
}


void Foam::mixingPlaneAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const scalarField& iF
) const
{
    fieldTransferBuffer_ = interfaceInternalField(iF);
}


Foam::tmp<Foam::labelField>
Foam::mixingPlaneAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList&
) const
{
    return shadowInterface().labelTransferBuffer();
}


Foam::tmp<Foam::scalarField>
Foam::mixingPlaneAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const scalarField&
) const
{
    return shadowInterface().fieldTransferBuffer();
}


// ************************************************************************* //
