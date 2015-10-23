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

\*---------------------------------------------------------------------------*/

#include "chemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::chemistryReader<ThermoType> >
Foam::chemistryReader<ThermoType>::New
(
    const dictionary& thermoDict
)
{
    // Let the chemistry reader type default to CHEMKIN
    // for backward compatability
    word chemistryReaderTypeName("chemkinReader");

    // otherwise use the specified reader
    thermoDict.readIfPresent("chemistryReader", chemistryReaderTypeName);

    Info<< "Selecting chemistryReader " << chemistryReaderTypeName << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(chemistryReaderTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "chemistryReader::New(const dictionary& thermoDict)"
        )   << "Unknown chemistryReader type "
            << chemistryReaderTypeName << nl << nl
            << "Valid  chemistryReaders are: " << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<chemistryReader<ThermoType> >(cstrIter()(thermoDict));
}


// ************************************************************************* //
