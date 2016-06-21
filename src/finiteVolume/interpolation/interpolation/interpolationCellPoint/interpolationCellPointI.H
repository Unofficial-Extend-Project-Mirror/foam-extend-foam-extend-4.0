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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::interpolationCellPoint<Type>::interpolate
(
    const cellPointWeight& cpw
) const
{
    const FixedList<scalar, 4>& weights = cpw.weights();
    const FixedList<label, 3>& faceVertices = cpw.faceVertices();

    Type t = psip_[faceVertices[0]]*weights[0];
    t += psip_[faceVertices[1]]*weights[1];
    t += psip_[faceVertices[2]]*weights[2];
    t += this->psi_[cpw.cell()]*weights[3];

    return t;
}


template<class Type>
inline Type Foam::interpolationCellPoint<Type>::interpolate
(
    const vector& position,
    const label celli,
    const label facei
) const
{
    return interpolate(cellPointWeight(this->pMesh_, position, celli, facei));
}


// ************************************************************************* //
