/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::IOField

Description
    A primitive field of type \<T\> with automated input and output.

SourceFiles
    IOField.C

\*---------------------------------------------------------------------------*/

#ifndef IOField_H
#define IOField_H

#include "regIOobject.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IOField
:
    public regIOobject,
    public Field<Type>
{

public:

    TypeName("Field");


    // Constructors

        //- Construct from IOobject
        IOField
        (
            const IOobject&
        );

        //- Construct from components
        IOField
        (
            const IOobject&,
            const Field<Type>&
        );

        //- Construct from IOobject and size (does not set values)
        IOField
        (
            const IOobject&,
            const label size
        );


    // Destructor

        virtual ~IOField();


    // Member functions

        bool writeData(Ostream& os) const
        {
            return (os << static_cast<const Field<Type>&>(*this)).good();
        }


    // Member operators

        void operator=(const IOField<Type>& iof)
        {
            Field<Type>::operator=(iof);
        }

        void operator=(const Field<Type>& f)
        {
            Field<Type>::operator=(f);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
