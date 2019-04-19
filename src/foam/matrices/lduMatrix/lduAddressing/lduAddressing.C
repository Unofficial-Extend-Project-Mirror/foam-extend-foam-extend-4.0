/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "lduAddressing.H"
#include "extendedLduAddressing.H"
#include "demandDrivenData.H"
#include "dynamicLabelList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::lduAddressing::calcLosort() const
{
    if (losortPtr_)
    {
        FatalErrorIn("lduAddressing::calcLosort() const")
            << "losort already calculated"
            << abort(FatalError);
    }

    // Scan the neighbour list to find out how many times the cell
    // appears as a neighbour of the face. Done this way to avoid guessing
    // and resizing list
    labelList nNbrOfFace(size(), 0);

    const unallocLabelList& nbr = upperAddr();

    forAll (nbr, nbrI)
    {
        nNbrOfFace[nbr[nbrI]]++;
    }

    // Create temporary neighbour addressing
    labelListList cellNbrFaces(size());

    forAll (cellNbrFaces, cellI)
    {
        cellNbrFaces[cellI].setSize(nNbrOfFace[cellI]);
    }

    // Reset the list of number of neighbours to zero
    nNbrOfFace = 0;

    // Scatter the neighbour faces
    forAll (nbr, nbrI)
    {
        cellNbrFaces[nbr[nbrI]][nNbrOfFace[nbr[nbrI]]] = nbrI;

        nNbrOfFace[nbr[nbrI]]++;
    }

    // Gather the neighbours into the losort array
    losortPtr_ = new labelList(nbr.size(), -1);

    labelList& lst = *losortPtr_;

    // Set counter for losort
    label lstI = 0;

    forAll (cellNbrFaces, cellI)
    {
        const labelList& curNbr = cellNbrFaces[cellI];

        forAll (curNbr, curNbrI)
        {
            lst[lstI] = curNbr[curNbrI];
            lstI++;
        }
    }
}


void Foam::lduAddressing::calcOwnerStart() const
{
    if (ownerStartPtr_)
    {
        FatalErrorIn("lduAddressing::calcOwnerStart() const")
            << "owner start already calculated"
            << abort(FatalError);
    }

    const labelList& own = lowerAddr();

    ownerStartPtr_ = new labelList(size() + 1, own.size());

    labelList& ownStart = *ownerStartPtr_;

    // Set up first lookup by hand
    ownStart[0] = 0;
    label nOwnStart = 0;
    label i = 1;

    forAll (own, faceI)
    {
        label curOwn = own[faceI];

        if (curOwn > nOwnStart)
        {
            while (i <= curOwn)
            {
                ownStart[i++] = faceI;
            }

            nOwnStart = curOwn;
        }
    }
}


void Foam::lduAddressing::calcLosortStart() const
{
    if (losortStartPtr_)
    {
        FatalErrorIn("lduAddressing::calcLosortStart() const")
            << "losort start already calculated"
            << abort(FatalError);
    }

    const labelList& nbr = upperAddr();

    losortStartPtr_ = new labelList(size() + 1, nbr.size());

    labelList& lsrtStart = *losortStartPtr_;

    const labelList& lsrt = losortAddr();

    // Set up first lookup by hand
    lsrtStart[0] = 0;
    label nLsrtStart = 0;
    label i = 0;

    forAll (lsrt, faceI)
    {
        // Get neighbour
        const label curNbr = nbr[lsrt[faceI]];

        if (curNbr > nLsrtStart)
        {
            while (i <= curNbr)
            {
                lsrtStart[i++] = faceI;
            }

            nLsrtStart = curNbr;
        }
    }

    // Set up last lookup by hand
    lsrtStart[size()] = nbr.size();
}


void Foam::lduAddressing::calcInternalBoundaryEqnCoeffs
(
    const lduInterfaceFieldPtrsList& lduInterfaces
) const
{
    if
    (
        internalEqnCoeffsPtr_
     || flippedInternalEqnCoeffsPtr_
     || boundaryEqnCoeffs_.size()
    )
    {
        FatalErrorIn("lduAddressing::calcInternalBoundaryCoeffs() const")
            << "Internal/boundary equation coefficients already calculated."
            << abort(FatalError);
    }

    // Get number of internal coefficients (number of faces)
    const label nInternalCoeffs = upperAddr().size();

    // Allocate insertion-friendly storage with enough memory
    dynamicLabelList internalCoeffs(nInternalCoeffs);
    dynamicLabelList flippedInternalCoeffs(nInternalCoeffs);

    // Initialise boundary equation coefficients for all coupled patches with
    // enough storage
    boundaryEqnCoeffs_.setSize(lduInterfaces.size());
    forAll (lduInterfaces, intI)
    {
        if (lduInterfaces.set(intI))
        {
            boundaryEqnCoeffs_.set
            (
                intI,
                new dynamicLabelList
                (
                    lduInterfaces[intI].coupledInterface().faceCells().size()
                )
            );
        }
    }

    // First, we need to mark boundary equations with associated interface
    // index. Note: this does not take into account corner cells that have
    // faces on more than one coupled interface. Don't care about them during
    // preconditioning at the moment. VV, 23/Jun/2017.
    labelList boundaryEqnInterfaceIndex(size_, -1);

    // Loop through interfaces
    forAll (lduInterfaces, intI)
    {
        // Check whether the interface is set
        if (lduInterfaces.set(intI))
        {
            // Get boundary equations/rows (face cells)
            const unallocLabelList& boundaryEqns =
                lduInterfaces[intI].coupledInterface().faceCells();

            // Loop through boundary equations and mark them
            forAll (boundaryEqns, beI)
            {
                boundaryEqnInterfaceIndex[boundaryEqns[beI]] = intI;
            }
        }
    }

    // Get lower/upper (owner/neighbour) addressing
    const unallocLabelList& own = lowerAddr();
    const unallocLabelList& nei = upperAddr();

    // Loop through upper triangle and filter coefficients (faces)
    forAll (own, coeffI)
    {
        // Get owner/neighbour (row/column) for this coefficient
        const label& ownI = own[coeffI];
        const label& neiI = nei[coeffI];

        if
        (
            boundaryEqnInterfaceIndex[ownI] != -1
         && boundaryEqnInterfaceIndex[neiI] != -1
        )
        {
            // Both owner and neigbour are at the boundary, append to boundary
            // coeffs list. Note: it is possible that owner/neighbour do not
            // have the same interface index since we ignored corner cells. Take
            // the owner index
            boundaryEqnCoeffs_[boundaryEqnInterfaceIndex[ownI]].append(coeffI);
        }
        else
        {
            // If owner is a boundary cell and neighbour is not a boundary cell,
            // we need to mark this face for flipping
            if
            (
                boundaryEqnInterfaceIndex[ownI] != -1
             && boundaryEqnInterfaceIndex[neiI] == -1
            )
            {
                flippedInternalCoeffs.append(coeffI);
            }
            else
            {
                internalCoeffs.append(coeffI);
            }
        }
    }

    // Reuse dynamic lists to initialise data members
    internalEqnCoeffsPtr_ = new labelList(internalCoeffs.xfer());
    flippedInternalEqnCoeffsPtr_ = new labelList(flippedInternalCoeffs.xfer());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduAddressing::lduAddressing(const label nEqns)
:
    size_(nEqns),
    losortPtr_(nullptr),
    ownerStartPtr_(nullptr),
    losortStartPtr_(nullptr),
    extendedAddr_(5),
    internalEqnCoeffsPtr_(nullptr),
    flippedInternalEqnCoeffsPtr_(nullptr),
    boundaryEqnCoeffs_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::lduAddressing::~lduAddressing()
{
    deleteDemandDrivenData(losortPtr_);
    deleteDemandDrivenData(ownerStartPtr_);
    deleteDemandDrivenData(losortStartPtr_);
    deleteDemandDrivenData(internalEqnCoeffsPtr_);
    deleteDemandDrivenData(flippedInternalEqnCoeffsPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::unallocLabelList& Foam::lduAddressing::losortAddr() const
{
    if (!losortPtr_)
    {
        calcLosort();
    }

    return *losortPtr_;
}


const Foam::unallocLabelList& Foam::lduAddressing::ownerStartAddr() const
{
    if (!ownerStartPtr_)
    {
        calcOwnerStart();
    }

    return *ownerStartPtr_;
}


const Foam::unallocLabelList& Foam::lduAddressing::losortStartAddr() const
{
    if (!losortStartPtr_)
    {
        calcLosortStart();
    }

    return *losortStartPtr_;
}


// Return edge index given owner and neighbour label
Foam::label Foam::lduAddressing::triIndex(const label a, const label b) const
{
    label own = min(a, b);

    label nbr = max(a, b);

    label startLabel = ownerStartAddr()[own];

    label endLabel = ownerStartAddr()[own + 1];

    const unallocLabelList& neighbour = upperAddr();

    for (label i = startLabel; i < endLabel; i++)
    {
        if (neighbour[i] == nbr)
        {
            return i;
        }
    }

    // If neighbour has not been found, something has gone seriously
    // wrong with the addressing mechanism
    FatalErrorIn
    (
        "lduAddressing::triIndex(const label owner, const label nbr) const"
    )   << "neighbour " << nbr << " not found for owner " << own << ". "
        << "Problem with addressing"
        << abort(FatalError);

    return -1;
}


const Foam::extendedLduAddressing&
Foam::lduAddressing::extendedAddr(const label p) const
{
    if (p == 0 || p > 4)
    {
        FatalErrorIn
        (
            "const Foam::extendedLduAddressing& "
            "Foam::lduAddressing::extendedAddr(const label p) const"
        )   << "Currently supported extended addressing fill-in only "
            << "between order 1 and 4"
            << abort(FatalError);
    }

    if (!extendedAddr_.set(p))
    {
        extendedAddr_.set(p, new extendedLduAddressing(*this, p));
    }

    return extendedAddr_[p];
}


const Foam::unallocLabelList& Foam::lduAddressing::internalEqnCoeffs
(
    const lduInterfaceFieldPtrsList& lduInterfaces
) const
{
    if (!internalEqnCoeffsPtr_)
    {
        calcInternalBoundaryEqnCoeffs(lduInterfaces);
    }

    return *internalEqnCoeffsPtr_;
}


const Foam::unallocLabelList& Foam::lduAddressing::flippedInternalEqnCoeffs
(
    const lduInterfaceFieldPtrsList& lduInterfaces
) const
{
    if (!flippedInternalEqnCoeffsPtr_)
    {
        calcInternalBoundaryEqnCoeffs(lduInterfaces);
    }

    return *flippedInternalEqnCoeffsPtr_;
}


const Foam::dynamicLabelList& Foam::lduAddressing::boundaryEqnCoeffs
(
    const lduInterfaceFieldPtrsList& lduInterfaces,
    const label intI
) const
{
    if (intI > lduInterfaces.size() - 1 || intI < 0)
    {
        FatalErrorIn
        (
            "const Foam::PtrList<labelList>& "
            "Foam::lduAddressing::boundaryEqnCoeffs"
            "\n("
            "\n    const lduInterfaceFieldPtrsList& lduInterfaces,"
            "\n    const label p"
            "\n) const"
        )   << "Invalid interface index specified: " << intI << nl
            << "Number of coupled interfaces: " << lduInterfaces.size()
            << abort(FatalError);
    }

    if (!boundaryEqnCoeffs_.set(intI))
    {
        calcInternalBoundaryEqnCoeffs(lduInterfaces);
    }

    return boundaryEqnCoeffs_[intI];
}


// ************************************************************************* //
