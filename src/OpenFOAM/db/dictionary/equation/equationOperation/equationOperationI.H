/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const equationOperation::sourceListType&
    equationOperation::sourceList() const
{
    return sourceList_;
}


inline equationOperation::sourceListType&
    equationOperation::sourceList()
{
    return sourceList_;
}


inline const label& equationOperation::sourceIndex() const
{
    return sourceIndex_;
}


inline label& equationOperation::sourceIndex()
{
    return sourceIndex_;
}

inline const label& equationOperation::dictLookupIndex() const
{
    return dictLookupIndex_;
}


inline label& equationOperation::dictLookupIndex()
{
    return dictLookupIndex_;
}


inline const equationOperation::operationType&
    equationOperation::operation() const
{
    return operation_;
}


inline equationOperation::operationType& equationOperation::operation()
{
    return operation_;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
