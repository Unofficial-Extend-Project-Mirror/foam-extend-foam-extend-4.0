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

#ifndef helperFunctionsFrontalMarking_H
#define helperFunctionsFrontalMarking_H

#include "labelList.H"
#include "LongList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace help
{

//- template implementation of the fron-marking algorithm
//- labelListType contains indices of elements satisfying the filtering criteria
//- startingIndex is the first element from where the search is started
//- neiOp determines which neighbouring elements are available. It requires
//- an operator neiOp(const label, DynList<label>&) to be defined
//- filterOp determines which neighbouring elements shall be stored
//- in the result and used in the front. It requires an operator
//- bool filterOp(const label) which return true if the element shall be part
//- of the front
template<class labelListType, class neiOp, class filterOp>
void frontalMarking
(
    labelListType& result,
    const label startingIndex,
    const neiOp& neighbourCalculator,
    const filterOp& selector
);

//- templated implementation of functionality for finding discrete groups
//- of elements depending on the user-defined criteria
//- the first argument is the result of the operation. Each element is assigned
//- to a group amd the ones failing the critearia are set to -1
//- neiOp is the functionality for determining neighbours of the currently
//- processed front element. This class requires an operator
//- neiOp(const label, DynList<label>&)
//- filterOp represents the user-defined criteria and shall have an operator
//- bool filterOp(const label)
template<class labelListType, class neiOp, class filterOp>
label groupMarking
(
    labelListType& elementInGroup,
    const neiOp& neighbourCalculator,
    const filterOp& selector
);

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "helperFunctionsFrontalMarking.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
