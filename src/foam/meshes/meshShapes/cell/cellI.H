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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
inline Foam::cell::cell()
{}


// Construct given size
inline Foam::cell::cell(label s)
:
    labelList(s, -1)
{}


// Construct from components
inline Foam::cell::cell(const UList<label>& lst)
:
    labelList(lst)
{}


inline Foam::cell::cell(const Xfer<labelList>& lst)
:
    labelList(lst)
{}


// Construct from Istream
inline Foam::cell::cell(Istream& is)
:
    labelList(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Number of faces
inline Foam::label Foam::cell::nFaces() const
{
    return size();
}


inline bool Foam::operator!=(const cell& a, const cell& b)
{
    return (!(a == b));
}

// ************************************************************************* //
