/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "FRWGraph.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, Foam::label width>
Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::FRWGraph<T, width>& DL
)
{
    os << DL.size() << "(" << endl;
    for(register label i=0;i<DL.size();++i)
    {
        os << width << "(";
        for(label j=0;j<width;++j)
            os << DL(i, j) << " " << endl;
        os << ")" << endl;
    }
    os << ")";
    return os;
}

/*
template<class T, Foam::label width>
Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::FRWGraph<T, width>& DL
)
{
    label size;
    T e;
    is >> size;
    DL.setSize(size);
    for(IndexType i=0;i<size;++i)
    {
        is >> e;
        DL[i] = e;
    }

    return is;
}
*/

// ************************************************************************* //
