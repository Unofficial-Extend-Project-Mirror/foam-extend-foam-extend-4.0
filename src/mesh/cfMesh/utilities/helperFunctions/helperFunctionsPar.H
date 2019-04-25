/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    helperFunctionsPar

Description
    Functions used for exchanging data between processors

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef helperFunctionsPar_H
#define helperFunctionsPar_H

#include "labelList.H"
#include "LongList.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace help
{

//- exchanges the data between processors such that in the end thwy all have
//- the necessary information to proceed with the local work
template<class sendOp, class combineOp, class T, class ListType>
void combineData(const sendOp&, combineOp&);

//- each processor exchanges the data with all of its neighbours
//- processors with lower labels first send the data to processors with
//- higher labels first and and then to the processors with lower labels
template<class T, class ListType, class scatterOp, class gatherOp>
void whisperReduce(const ListType&, const scatterOp&, gatherOp&);

//- send the data stored in the map to other processors and collects the data
//- sent from other processors into the list
template<class T, class ListType>
void exchangeMap
(
    const std::map<label, ListType>&,
    LongList<T>&,
    const Pstream::commsTypes commsType = Pstream::blocking
);

//- sends the data stored in a map to other processors and receives the data
//- received from other processors into another map
template<class T, class ListType>
void exchangeMap(const std::map<label, ListType>&, std::map<label, List<T> >&);

//- calculates the reverse addressing of the graph by transposing the graph
template<class RowType, template<class ListTypeArg> class GraphType>
void reverseAddressingSMP(const GraphType<RowType>&, GraphType<RowType>&);

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "helperFunctionsPar.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
