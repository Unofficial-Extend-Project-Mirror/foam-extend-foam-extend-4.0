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

Class
    Foam::IOPtrList

Description
    A PtrList of objects of type \<T\> with automated input and output.

SourceFiles
    IOPtrList.C

\*---------------------------------------------------------------------------*/

#ifndef IOPtrList_H
#define IOPtrList_H

#include "PtrList.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOPtrList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class IOPtrList
:
    public regIOobject,
    public PtrList<T>
{

public:

    //- Runtime type information
    TypeName("PtrList");


    // Constructors

        //- Construct from IOobject using given Istream constructor class
        template<class INew>
        IOPtrList(const IOobject&, const INew&);

        //- Construct from IOobject
        IOPtrList(const IOobject&);

        //- Construct from IOobject with given size
        IOPtrList(const IOobject&, const label);

        //- Construct from IOobject and a PtrList
        IOPtrList(const IOobject&, const PtrList<T>&);

        //- Construct by transferring the PtrList contents
        IOPtrList(const IOobject&, const Xfer<PtrList<T> >&);


    // Destructor

        virtual ~IOPtrList();


    // Member functions

        bool writeData(Ostream&) const;


    // Member operators

        void operator=(const IOPtrList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOPtrList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
