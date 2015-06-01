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

Typedef
    Foam::tuple2List

Description
    List of paired values, with sorting capability

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef tuple2Lists_H
#define tuple2Lists_H

#include "Tuple2.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef List<Tuple2<label, label> > labelLabelList;
    typedef List<Tuple2<label, scalar> > labelScalarList;
    typedef List<Tuple2<scalar, label> > scalarLabelList;
    typedef List<Tuple2<scalar, scalar> > scalarScalarList;

template<class type1, class type2>
bool lessFirst(Tuple2<type1, type2> t1, Tuple2<type1, type2> t2)
{
    return (t1.first() < t2.first());
}

template<class type1, class type2>
bool lessFirstSecond(Tuple2<type1, type2> t1, Tuple2<type1, type2> t2)
{
    return
    (
        (t1.first() < t2.first())
     || ((t1.first() == t2.first()) && (t1.second() < t2.second()))
    );
}

template<class type1, class type2>
bool lessSecond(Tuple2<type1, type2> t1, Tuple2<type1, type2> t2)
{
    return (t1.second() < t2.second());
}

template<class type1, class type2>
bool lessSecondFirst(Tuple2<type1, type2> t1, Tuple2<type1, type2> t2)
{
    return
    (
        (t2.first() < t1.first())
     || ((t2.first() == t1.first()) && (t2.second() < t1.second()))
    );
}

template<class type1, class type2>
void sortTuple2ListBy1st(List<Tuple2<type1, type2> >& t2l)
{
    sort(t2l, lessFirst<type1, type2>);
}

template<class type1, class type2>
void sortTuple2ListBy2nd(List<Tuple2<type1, type2> >& t2l)
{
    sort(t2l, lessSecond<type1, type2>);
}

template<class type1, class type2>
void sortTuple2ListBy1stThen2nd(List<Tuple2<type1, type2> >& t2l)
{
    sort(t2l, lessFirstSecond<type1, type2>);
}

template<class type1, class type2>
void sortTuple2ListBy2ndThen1st(List<Tuple2<type1, type2> >& t2l)
{
    sort(t2l, lessSecondFirst<type1, type2>);
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
