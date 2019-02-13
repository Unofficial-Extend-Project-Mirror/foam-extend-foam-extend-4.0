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
    faceListPMG

Description
    This is a container with additional size to prevent re-allocation
    every time it is resized

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef faceListPMG_H
#define faceListPMG_H

#include "regIOobject.H"
#include "faceList.H"
#include "Istream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class faceListPMG
:   public regIOobject,
    public faceList
{
    // Private data
        //- number of used elements
        label nElmts_;

    // Disallow bitwise assignment
        void operator=(const faceListPMG&);

        faceListPMG(const faceListPMG&);

    // Disallow transfer from faceList
        void transfer(faceList&);

public:

        TypeName("faceList");

    // Constructors
        //- construct from IOobject
        inline faceListPMG(const IOobject&);

        //- construct from IOobject and size
        inline faceListPMG(const IOobject&, const label);

        //- construct from IOobject and faceList
        inline faceListPMG(const IOobject&, const faceList&);

    // Destructor
        inline ~faceListPMG();

    // Member functions
        //- return the number of used elements
        inline label size() const;

        //- set the number of used elements
        inline void setSize(const label nElmts);

        //- set the size to zero
        inline void clear();

        //- add a face at the end of the list
        inline void append(const face&);

        //- return an element with bound checking
        inline face& newElmt(const label);

        //- read/write the list onto disk
        inline bool writeData(Ostream&) const;

    // Member operators
        inline void operator=(const faceList&);

        friend inline Ostream& operator<<(Ostream&, const faceListPMG&);

        friend inline Istream& operator>>(Istream&, faceListPMG&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceListPMGI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
