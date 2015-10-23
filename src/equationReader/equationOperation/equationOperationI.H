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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

inline const equationOperation::sourceTypeEnum&
    equationOperation::sourceType() const
{
    return source_;
}


inline equationOperation::sourceTypeEnum&
    equationOperation::sourceType()
{
    return source_;
}


inline const label& equationOperation::sourceIndex() const
{
    return sourceIndex_;
}


inline label& equationOperation::sourceIndex()
{
    return sourceIndex_;
}


inline const label& equationOperation::componentIndex() const
{
    return componentIndex_;
}


inline label& equationOperation::componentIndex()
{
    return componentIndex_;
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
