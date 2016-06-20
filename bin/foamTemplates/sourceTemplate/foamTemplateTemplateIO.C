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

#include "ClassName.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::ClassName<TemplateArgument>::ClassName(Istream& is)
:
    base1(is),
    base2(is),
    member1(is),
    member2(is)
{
    // Check state of Istream
    is.check("Foam::ClassName<TemplateArgument>::ClassName(Foam::Istream&)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    ClassName<TemplateArgument>&
)
{
    // Check state of Istream
    is.check
    (
        "Foam::Istream& Foam::operator>>"
        "(Foam::Istream&, Foam::ClassName<TemplateArgument>&)"
    );

    return is;
}


template<TemplateClassArgument>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ClassName<TemplateArgument>&
)
{
    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Ostream&, const ClassName<TemplateArgument>&)"
    );

    return os;
}


// ************************************************************************* //
