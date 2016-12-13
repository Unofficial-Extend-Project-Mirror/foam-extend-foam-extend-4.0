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

#include "extendedLduAddressing.H"
#include "mapPolyMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(extendedLduAddressing, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::extendedLduAddressing::clearOut() const
{
    deleteDemandDrivenData(extendedLowerPtr_);
    deleteDemandDrivenData(extendedUpperPtr_);
    deleteDemandDrivenData(faceMapPtr_);
    deleteDemandDrivenData(extendedLosortPtr_);
    deleteDemandDrivenData(extendedOwnerStartPtr_);
    deleteDemandDrivenData(extendedLosortStartPtr_);
}


void Foam::extendedLduAddressing::markNeighbours
(
    const label& masterCellI,
    const labelList& nbrCells,
    const labelListList& cellCells,
    HashSet<label>& extendedNbrs,
    label curLevel
) const
{
    if (curLevel > -1)
    {
        // For each neighbouring cell, loop through their neighbours
        forAll (nbrCells, i)
        {
            // Get neighbouring cell index
            const label& nbrCellI = nbrCells[i];

            // Insert this neighbouring cell in the extended cell neighbours
            // hash set if its index is greater than current cell
            if (nbrCellI > masterCellI)
            {
                extendedNbrs.insert(nbrCellI);
            }

            // Recursive call with decremented level and neighbours of this
            // neighbouring cells
            markNeighbours
            (
                masterCellI,
                cellCells[nbrCellI],
                cellCells,
                extendedNbrs,
                curLevel - 1
            );
        }
    }
}


void Foam::extendedLduAddressing::calcCellCells
(
    labelListList& cellCellAddr
) const
{
    // Algorithm taken from primiteMeshCellCells.C in calcCellCells() member
    // function.
    labelList ncc(lduAddr_.size(), 0);

    const unallocLabelList& own = lowerAddr();
    const unallocLabelList& nei = upperAddr();

    forAll (nei, faceI)
    {
        ncc[own[faceI]]++;
        ncc[nei[faceI]]++;
    }

    // Size the list
    cellCellAddr.setSize(ncc.size());

    forAll (cellCellAddr, cellI)
    {
        cellCellAddr[cellI].setSize(ncc[cellI]);
    }
    ncc = 0;

    forAll (nei, faceI)
    {
        label ownCellI = own[faceI];
        label neiCellI = nei[faceI];

        cellCellAddr[ownCellI][ncc[ownCellI]++] = neiCellI;
        cellCellAddr[neiCellI][ncc[neiCellI]++] = ownCellI;
    }
}


void Foam::extendedLduAddressing::calcExtendedLowerUpper() const
{
    if (extendedLowerPtr_ && extendedUpperPtr_)
    {
        FatalErrorIn
        (
            "extendedLduAddressing::calcExtendeLowerUpper() const"
        )
            << "Extended lower/upper addressing are already calculated"
            << abort(FatalError);
    }
    else if (extendedLowerPtr_ || extendedUpperPtr_)
    {
        FatalErrorIn
        (
            "extendedLduAddressing::calcExtendeLowerUpper() const"
        )
            << "Extended lower/upper addressing partially calculated."
            << "Something went terribly wrong. "
            << abort(FatalError);
    }

    // Allocate dynamic lists for extended owner/neighbour, which are later
    // used to define ordinary labelLists (extendedLower, extendedUpper)
    // Helper type definition
    typedef dynamicLabelList DynamicLabelList;

    // Get the number of faces
    const label nFaces = upperAddr().size();

    // Allocate extended owner/neighbour with 4*nFaces*extensionLevel as a
    // guess (based on hex mesh assumption) to prevent excessive resizing
    DynamicLabelList dExtOwn(4*nFaces*p_);
    DynamicLabelList dExtNei(4*nFaces*p_);

    // Create working hash set of extended neighbours for a given cell
    HashSet<label> extCellNbrs(32*p_);

    // Get number of cells in the mesh
    const label nCells = lduAddr_.size();

    // Get a list of cells neighbouring each cell.
    // Note optimisation to avoid a copy
    labelListList cellCells;
    calcCellCells(cellCells);

    // Loop through cells
    for (label cellI = 0; cellI < nCells; ++cellI)
    {
        // Get neighbouring cells of this cell
        const labelList& nbrCells = cellCells[cellI];

        // Mark the neighbours based on extension level
        markNeighbours(cellI, nbrCells, cellCells, extCellNbrs, p_);

        // Get sorted list of extended neighbours from the hash set. Warning:
        // .sortedToc() function returns List<key> by value
        const labelList sortedExtNbrs = extCellNbrs.sortedToc();

        // Loop through sorted extended neighbours
        forAll (sortedExtNbrs, extNbrI)
        {
            // Append cellI in the owner list and extNbrI in neighbour list
            dExtOwn.append(cellI);
            dExtNei.append(sortedExtNbrs[extNbrI]);
        }

        // Clear the hash set, keeping the allocated memory
        extCellNbrs.clear();
    }

    // Create extended lower/upper addressing, re-using dynamic lists
    extendedLowerPtr_ = new labelList(dExtOwn, true);
    extendedUpperPtr_ = new labelList(dExtNei, true);
}


void Foam::extendedLduAddressing::calcFaceMap() const
{
    if (faceMapPtr_)
    {
        FatalErrorIn("extendedLduAddressing::calcFaceMap() const")
            << "faceMap already calculated"
            << abort(FatalError);
    }

    // Get reference to ordinary owner/neighbour addressing
    const unallocLabelList& own = lowerAddr();
    const unallocLabelList& nbr = upperAddr();

    // Allocate memory for faceMap
    faceMapPtr_ = new labelList(own.size(), -1);
    labelList& faceMap = *faceMapPtr_;

    // Get reference to extended neighbour and owner start addressing
    const unallocLabelList& extNbr = extendedUpperAddr();
    const unallocLabelList& extOwnStart = extendedOwnerStartAddr();

    // Loop through ordinary faces
    forAll (nbr, faceI)
    {
        // Get owner and neighbour of this face
        const label& ownI = own[faceI];
        const label& nbrI = nbr[faceI];

        // Get start and end labels for this owner cell
        const label& startLabel = extOwnStart[ownI];
        const label& endLabel = extOwnStart[ownI + 1];

        // Loop through extended neighbouring faces (upper matrix coeffs)
        for
        (
            register label extFaceI = startLabel;
            extFaceI < endLabel;
            ++extFaceI
        )
        {
            if (extNbr[extFaceI] == nbrI)
            {
                faceMap[faceI] = extFaceI;
            }
        }
    }
}


void Foam::extendedLduAddressing::calcExtendedLosort() const
{
    if (extendedLosortPtr_)
    {
        FatalErrorIn("extendedLduAddressing::calcExtendedLosort() const")
            << "extendedLosort already calculated"
            << abort(FatalError);
    }

    // Algorithm is the same as the one in lduAddressing, the only difference
    // being usage of extended owner and neighbour instead of owner/neighbour
    // lists

    // Scan the extended neighbour list to find out how many times the cell
    // appears as a neighbour of the face. Done this way to avoid guessing
    // and resizing list
    const label matrixSize = lduAddr_.size();
    labelList nNbrOfFace(matrixSize, 0);

    const unallocLabelList& nbr = extendedUpperAddr();

    forAll (nbr, nbrI)
    {
        ++nNbrOfFace[nbr[nbrI]];
    }

    // Create temporary neighbour addressing
    labelListList cellNbrFaces(matrixSize);

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

        ++nNbrOfFace[nbr[nbrI]];
    }

    // Gather the neighbours into the losort array
    extendedLosortPtr_ = new labelList(nbr.size(), -1);

    labelList& lst = *extendedLosortPtr_;

    // Set counter for losort
    label lstI = 0;

    forAll (cellNbrFaces, cellI)
    {
        const labelList& curNbr = cellNbrFaces[cellI];

        forAll (curNbr, curNbrI)
        {
            lst[lstI] = curNbr[curNbrI];
            ++lstI;
        }
    }
}


void Foam::extendedLduAddressing::calcExtendedOwnerStart() const
{
    if (extendedOwnerStartPtr_)
    {
        FatalErrorIn
        (
            "extendedLduAddressing::calcExtendedOwnerStart() const"
        )
            << "extended owner start already calculated"
            << abort(FatalError);
    }

    // Algorithm is the same as the one in lduAddressing, the only difference
    // being usage of extended owner and neighbour instead of owner/neighbour
    // lists

    const labelList& own = extendedLowerAddr();

    extendedOwnerStartPtr_ = new labelList(lduAddr().size() + 1, own.size());

    labelList& ownStart = *extendedOwnerStartPtr_;

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


void Foam::extendedLduAddressing::calcExtendedLosortStart() const
{
    if (extendedLosortStartPtr_)
    {
        FatalErrorIn
        (
            "extendedLduAddressing::calcExtendedLosortStart() const"
        )
            << "extended losort start already calculated"
            << abort(FatalError);
    }

    // Algorithm is the same as the one in lduAddressing, the only difference
    // being usage of extended owner and neighbour instead of owner/neighbour
    // lists

    extendedLosortStartPtr_ = new labelList(lduAddr().size() + 1, 0);

    labelList& lsrtStart = *extendedLosortStartPtr_;

    const labelList& nbr = extendedUpperAddr();

    const labelList& lsrt = extendedLosortAddr();

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
    lsrtStart[lduAddr_.size()] = nbr.size();
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::extendedLduAddressing::extendedLduAddressing
(
    const lduAddressing& lduAddr,
    const label extensionLevel
)
:
    lduAddr_(lduAddr),
    p_(extensionLevel),
    extendedLowerPtr_(NULL),
    extendedUpperPtr_(NULL),
    faceMapPtr_(NULL),
    extendedLosortPtr_(NULL),
    extendedOwnerStartPtr_(NULL),
    extendedLosortStartPtr_(NULL)
{
    // Issue an error if a negative extension level is selected
    if (p_ < 0)
    {
        FatalErrorIn
        (
            "extendedLduAddressing::extendedLduAddressing"
        )
            << "Negative extension level not allowed."
            << abort(FatalError);
    }
    // Disallow extension level 0 as it is the same as ordinary lduAddressing
    else if (p_ == 0)
    {
        FatalErrorIn
        (
            "extendedLduAddressing::extendedLduAddressing"
        )
            << "Extension level 0 not allowed as it is the same as ordinary "
            << "lduAddressing."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * /

Foam::extendedLduAddressing::~extendedLduAddressing()
{
    clearOut();
}


// * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

const Foam::unallocLabelList&
Foam::extendedLduAddressing::extendedLowerAddr() const
{
    // Extended lower and upper are calculated together. If one of
    // these lists has been calculated, and one of them hasn't been
    // calculated, something went terribly wrong with addressing
    if(!extendedLowerPtr_ || !extendedUpperPtr_)
    {
        calcExtendedLowerUpper();
    }

    return *extendedLowerPtr_;
}


const Foam::unallocLabelList&
Foam::extendedLduAddressing::extendedUpperAddr() const
{
    // Extended lower and upper are calculated together. If one of
    // these lists has been calculated, and one of them hasn't been
    // calculated, something went terribly wrong with addressing
    if(!extendedLowerPtr_ || !extendedUpperPtr_)
    {
        calcExtendedLowerUpper();
    }

    return *extendedUpperPtr_;
}


const Foam::unallocLabelList& Foam::extendedLduAddressing::faceMap() const
{
    if(!faceMapPtr_)
    {
        calcFaceMap();
    }

    return *faceMapPtr_;
}


const Foam::unallocLabelList&
Foam::extendedLduAddressing::extendedLosortAddr() const
{
    if (!extendedLosortPtr_)
    {
        calcExtendedLosort();
    }

    return *extendedLosortPtr_;
}


const Foam::unallocLabelList&
Foam::extendedLduAddressing::extendedOwnerStartAddr() const
{
    if (!extendedOwnerStartPtr_)
    {
        calcExtendedOwnerStart();
    }

    return *extendedOwnerStartPtr_;
}


const Foam::unallocLabelList&
Foam::extendedLduAddressing::extendedLosortStartAddr() const
{
    if (!extendedLosortStartPtr_)
    {
        calcExtendedLosortStart();
    }

    return *extendedLosortStartPtr_;
}


Foam::label
Foam::extendedLduAddressing::extendedTriIndex
(
    const label a,
    const label b
) const
{
    label own = min(a, b);

    label nbr = max(a, b);

    label startLabel = extendedOwnerStartAddr()[own];

    label endLabel = extendedOwnerStartAddr()[own + 1];

    const unallocLabelList& neighbour = extendedUpperAddr();

    for (label i = startLabel; i < endLabel; ++i)
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
        "extendedLduAddressing::extendedTriIndex\n"
        "(\n"
        "   const label owner,\n"
        "   const label neighbour\n"
        ") const \n"
    )   << "neighbour " << nbr << " not found for owner " << own << ". "
        << "Problem with extended addressing"
        << abort(FatalError);

    return -1;
}


bool Foam::extendedLduAddressing::updateMesh(const mapPolyMesh& mpm) const
{
    if (debug)
    {
        InfoIn
        (
            "extendedLduAddressing::updateMesh(const mapPolyMesh&) const"
        )   << "Clearing extendedLduAddressing data" << endl;
    }

    clearOut();

    return true;
}


// ************************************************************************* //
