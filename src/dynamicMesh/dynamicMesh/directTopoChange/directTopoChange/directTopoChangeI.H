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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline bool Foam::directTopoChange::pointRemoved(const label pointI) const
{
    const point& pt = points_[pointI];

    return
        pt.x() > 0.5*greatPoint.x()
     && pt.y() > 0.5*greatPoint.y()
     && pt.z() > 0.5*greatPoint.z();
}


inline bool Foam::directTopoChange::faceRemoved(const label faceI) const
{
    return faces_[faceI].size() == 0;
}


inline bool Foam::directTopoChange::cellRemoved(const label cellI) const
{
    return cellMap_[cellI] == -2;
}


inline void Foam::directTopoChange::setNumPatches(const label nPatches)
{
    nPatches_ = nPatches;
}


// ************************************************************************* //
