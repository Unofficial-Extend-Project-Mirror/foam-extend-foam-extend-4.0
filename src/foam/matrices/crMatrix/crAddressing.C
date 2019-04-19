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

Class
    crAddressing

Description
    Compressed row addressing class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "crAddressing.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::crAddressing::setRowSizes(const labelList& rowSizes)
{
    if (rowSizes.size() != nRows_)
    {
        FatalErrorIn
        (
            "void crAddressing::setRowSizes(const labelList& rowSizes)"
        )   << "Incorrect size of rowSizes: nRows =" << nRows_
            << " rowSizes = " << rowSizes.size()
            << abort(FatalError);

    }

    // Accumulate rowSizes
    rowStart_[0] = 0;

    forAll (rowSizes, i)
    {
        rowStart_[i + 1] = rowStart_[i] + rowSizes[i];
    }

    // Resize and clear column array
    column_.setSize(rowStart_[nRows_]);
    column_ = 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
Foam::crAddressing::crAddressing()
:
    refCount(),
    nRows_(0),
    nCols_(0),
    rowStart_(1),
    column_(0)
{
    rowStart_[0] = 1;
}


// Construct from given size
Foam::crAddressing::crAddressing
(
    const label nRows,
    const label nCols,
    const labelList& rowSizes
)
:
    refCount(),
    nRows_(nRows),
    nCols_(nCols),
    rowStart_(nRows + 1),
    column_(0)
{
    setRowSizes(rowSizes);
}


// Construct from components
Foam::crAddressing::crAddressing
(
    const label nRows,
    const label nCols,
    const labelList& row,
    const labelList& col
)
:
    refCount(),
    nRows_(nRows),
    nCols_(nCols),
    rowStart_(row),
    column_(col)
{}


// Construct as copy
Foam::crAddressing::crAddressing(const crAddressing& a)
:
    refCount(),
    nRows_(a.nRows_),
    nCols_(a.nCols_),
    rowStart_(a.rowStart_),
    column_(a.column_)
{}


// Construct from Istream
Foam::crAddressing::crAddressing(Istream& is)
:
    refCount(),
    nRows_(readLabel(is)),
    nCols_(readLabel(is)),
    rowStart_(is),
    column_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::crAddressing> Foam::crAddressing::T() const
{
    const labelList& myRow = rowStart();
    const labelList& myCol = column();

    labelList trRowSizes(nCols(), 0);

    // Count number of entries in a row
    for (label i = 0; i < nRows(); i++)
    {
        for (label ip = myRow[i]; ip < myRow[i + 1]; ip++)
        {
            trRowSizes[myCol[ip]]++;
        }
    }

    // Create transpose addressing
    tmp<crAddressing> ttranspose
    (
        new crAddressing(nCols(), nRows(), trRowSizes)
    );
    crAddressing& transpose = ttranspose();

    // Set coefficients
    const labelList& trRow = transpose.rowStart();
    labelList& trCol = transpose.column();
    trCol = 0;

    // Reset rowSizes to use as counter
    trRowSizes = 0;

    label j;

    for (label i = 0; i < nRows(); i++)
    {
        for (label ip = myRow[i]; ip < myRow[i + 1]; ip++)
        {
            j = myCol[ip];

            trCol[trRow[j] + trRowSizes[j]] = i;

            trRowSizes[j]++;
        }
    }

    return ttranspose;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::crAddressing::operator=(const crAddressing& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::crAddressing::operator=(const Foam::crAddressing&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    nRows_ = rhs.nRows_;
    nCols_ = rhs.nCols_;
    rowStart_ = rhs.rowStart_;
    column_ = rhs.column_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const crAddressing& a)
{
    os  << a.nRows_ << tab << a.nCols_ << nl
        << a.rowStart_ << nl
        << a.column_ << endl;

    return os;
}


// ************************************************************************* //
