/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "TableFile.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::TableFile<Type>::TableFile
(
    const word& entryName,
    const dictionary& dict
)
:
    TableBase<Type>(entryName, dict.subDict(entryName + "Coeffs")),
    fName_("none")
{
    const dictionary coeffs(dict.subDict(entryName + "Coeffs"));
    coeffs.lookup("fileName") >> fName_;

    fileName expandedFile(fName_);
    IFstream is(expandedFile.expand());

    if (!is.good())
    {
        FatalIOErrorInFunction
        (
            is
        )   << "Cannot open file." << exit(FatalIOError);
    }

    is  >> this->table_;

    TableBase<Type>::check();
}


template<class Type>
Foam::Function1Types::TableFile<Type>::TableFile(const TableFile<Type>& tbl)
:
    TableBase<Type>(tbl),
    fName_(tbl.fName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::TableFile<Type>::~TableFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::TableFile<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::END_STATEMENT << nl
        << indent << word(this->name() + "Coeffs") << nl
        << indent << token::BEGIN_BLOCK << nl << incrIndent;

    // Note: for TableBase write the dictionary entries it needs but not
    // the values themselves
    TableBase<Type>::writeEntries(os);

    os.writeKeyword("fileName")<< fName_ << token::END_STATEMENT << nl;
    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
