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

#include "fieldValue.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::fieldValue::name() const
{
    return name_;
}


inline const Foam::objectRegistry& Foam::fieldValue::obr() const
{
    return obr_;
}


inline bool Foam::fieldValue::active() const
{
    return active_;
}


inline const Foam::Switch& Foam::fieldValue::log() const
{
    return log_;
}


inline const Foam::word& Foam::fieldValue::sourceName() const
{
    return sourceName_;
}


inline const Foam::wordList& Foam::fieldValue::fields() const
{
    return fields_;
}


inline const Foam::Switch& Foam::fieldValue::valueOutput() const
{
    return valueOutput_;
}


inline const Foam::fvMesh& Foam::fieldValue::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// ************************************************************************* //
