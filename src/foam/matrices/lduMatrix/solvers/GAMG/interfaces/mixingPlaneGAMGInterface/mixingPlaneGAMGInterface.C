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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Contributor
    Martin Beaudoin, Hydro-Quebec, 2009.

\*---------------------------------------------------------------------------*/

#include "mixingPlaneGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingPlaneGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        mixingPlaneGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingPlaneGAMGInterface::mixingPlaneGAMGInterface
(
    const lduPrimitiveMesh& lduMesh,
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
:
    GAMGInterface(lduMesh),
    fineMixingPlaneInterface_
    (
        refCast<const mixingPlaneLduInterface>(fineInterface)
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::mixingPlaneGAMGInterface::~mixingPlaneGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::mixingPlaneGAMGInterface::agglomerateCoeffs
(
    const scalarField& fineCoeffs
) const
{
    // AMG agglomeration missing
    notImplemented("mixingPlaneGAMGInterface::agglomerateCoeffs");
    tmp<scalarField> tcoarseCoeffs(new scalarField(size(), 0.0));

    return tcoarseCoeffs;
}


bool Foam::mixingPlaneGAMGInterface::master() const
{
    return fineMixingPlaneInterface_.master();
}


Foam::label Foam::mixingPlaneGAMGInterface::shadowIndex() const
{
    return fineMixingPlaneInterface_.shadowIndex();
}


const Foam::mixingPlaneLduInterface&
Foam::mixingPlaneGAMGInterface::shadowInterface() const
{
    return refCast<const mixingPlaneLduInterface>
    (
        ldu().interfaces()[shadowIndex()]
    );
}


const Foam::labelListList& Foam::mixingPlaneGAMGInterface::addressing() const
{
    FatalErrorIn
    (
        "const labelListList& mixingPlaneGAMGInterface::addressing() const"
    )   << "Requested fine addressing at coarse level"
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::mixingPlaneGAMGInterface::weights() const
{
    FatalErrorIn
    (
        "const labelListList& mixingPlaneGAMGInterface::weights() const"
    )   << "Requested fine addressing at coarse level"
        << abort(FatalError);
    return scalarListList::null();
}


const Foam::tensorField& Foam::mixingPlaneGAMGInterface::forwardT() const
{
    return fineMixingPlaneInterface_.forwardT();
}


const Foam::tensorField& Foam::mixingPlaneGAMGInterface::reverseT() const
{
    return fineMixingPlaneInterface_.reverseT();
}


void Foam::mixingPlaneGAMGInterface::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    labelTransferBuffer_ = interfaceData;
}


Foam::tmp<Foam::labelField> Foam::mixingPlaneGAMGInterface::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    return this->shadowInterface().labelTransferBuffer();
}


void Foam::mixingPlaneGAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    labelTransferBuffer_ = interfaceInternalField(iF);
}


void Foam::mixingPlaneGAMGInterface::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const scalarField& iF
) const
{
    fieldTransferBuffer_ = interfaceInternalField(iF);
}


Foam::tmp<Foam::labelField>
Foam::mixingPlaneGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList&
) const
{
    return shadowInterface().labelTransferBuffer();
}


Foam::tmp<Foam::scalarField>
Foam::mixingPlaneGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const scalarField&
) const
{
    return shadowInterface().fieldTransferBuffer();
}


// ************************************************************************* //
