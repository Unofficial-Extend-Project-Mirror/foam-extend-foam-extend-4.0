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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    crAddressing

Description
    Compressed row addressing class

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "crAddressing.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::crAddressing::setRowCount(const labelList& count)
{
    if (count.size() != nRows_)
    {
        FatalErrorIn("void crAddressing::setRowCount(const labelList& count)")
            << "Incorrect size of count: nRows =" << nRows_
            << " count = " << count.size()
            << abort(FatalError);

    }

    // Accumulate count
    row_[0] = 0;

    forAll (count, i)
    {
        row_[i + 1] = row_[i] + count[i];
    }

    // Resize and clear column array
    col_.setSize(row_[nRows_]);
    col_ = 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from given size
Foam::crAddressing::crAddressing
(
    const label nRows,
    const label nCols,
    const labelList& count
)
:
    refCount(),
    nRows_(nRows),
    nCols_(nCols),
    row_(nRows + 1),
    col_(0)
{
    setRowCount(count);
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
    row_(row),
    col_(col)
{}


// Construct as copy
Foam::crAddressing::crAddressing(const crAddressing& a)
:
    refCount(),
    nRows_(a.nRows_),
    nCols_(a.nCols_),
    row_(a.row_),
    col_(a.col_)
{}


// Construct from Istream
Foam::crAddressing::crAddressing(Istream& is)
:
    refCount(),
    nRows_(readLabel(is)),
    nCols_(readLabel(is)),
    row_(is),
    col_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::crAddressing> Foam::crAddressing::T() const
{
    const labelList& myRow = row();
    const labelList& myCol = col();

    labelList trCount(nCols(), 0);

    // Count number of entries in a row
    for (label i = 0; i < nRows(); i++)
    {
        for (label ip = myRow[i]; ip < myRow[i + 1]; ip++)
        {
            trCount[myCol[ip]]++;
        }
    }

    // Create transpose addressing
    tmp<crAddressing> ttranspose(new crAddressing(nCols(), nRows(), trCount));
    crAddressing& transpose = ttranspose();

    // Set coefficients
    const labelList& trRow = transpose.row();
    labelList& trCol = transpose.col();
    trCol = 0;

    // Reset count to use as counter
    trCount = 0;

    label j;

    for (label i = 0; i < nRows(); i++)
    {
        for (label ip = myRow[i]; ip < myRow[i + 1]; ip++)
        {
            j = myCol[ip];

            trCol[trRow[j] + trCount[j]] = i;

            trCount[j]++;
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
        FatalErrorIn("Foam::crAddressing::operator=(const Foam::crAddressing&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    nRows_ = rhs.nRows_;
    nCols_ = rhs.nCols_;
    row_ = rhs.row_;
    col_ = rhs.col_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const crAddressing& a)
{
    os  << a.nRows_ << tab << a.nCols_ << nl
        << a.row_ << nl
        << a.col_ << endl;

    return os;
}


// ************************************************************************* //
