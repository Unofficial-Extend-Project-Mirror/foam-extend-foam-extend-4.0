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

Class
    Foam::timeClusterList

Description
    List of timeClusters with some additional functions

SourceFiles
    timeClusterList.C

\*---------------------------------------------------------------------------*/

#ifndef timeClusterList_H
#define timeClusterList_H

#include "timeCluster.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class timeClusterList Declaration
\*---------------------------------------------------------------------------*/

class timeClusterList
:
    public List<timeCluster>
{

public:

    // Static data members

        static const char* const typeName;

    // Constructors

        //- Null constructor
        timeClusterList();

        //- Construct with given size.
        explicit timeClusterList(const label);

        //- Construct from a timeCluster.
        timeClusterList(const timeCluster&);

        //- Construct with given size and value for all elements.
        timeClusterList(const label, const timeCluster&);

        //- Construct as a sublist of another timeClusterList. This is not
        //  like the SubList in that this is not a new class, rather it is only
        //  producing a new copy that contains the specified entries.
        timeClusterList(const labelList&, const timeClusterList&);

        //- Construct from Istream.
        timeClusterList(Istream&);

    // Member functions

        //- Sort by global time
        void globalSort();

        //- Append another timeClusterList on the end of this one
        void append(const timeClusterList& tclIn);

        //- Append a timeCluster on the end of this list
        void append(const timeCluster& tcIn);

        //- Remove timeClusters with empty instantLists return false if all are
        //  empty
        bool purgeEmpties();

        //- Return a sublist using the specified indexes
        timeClusterList selectiveSubList(const labelList&) const;

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
