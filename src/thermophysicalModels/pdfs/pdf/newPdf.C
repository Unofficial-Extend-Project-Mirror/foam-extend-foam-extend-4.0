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

\*---------------------------------------------------------------------------*/

#include "pdf.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

namespace Foam 
{

autoPtr<Foam::pdf> Foam::pdf::New
(
    const dictionary& dict,
    Random& rndGen
)
{
    word pdfType
    (
        dict.lookup("pdfType")
    );

        Info << "Selecting pdfType "
            << pdfType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(pdfType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "pdf::New(const dictionary&, Random&) : " << endl
            << "    unknown pdfType type "
            << pdfType
            << ", constructor not in hash table" << endl << endl
            << "    Valid pdf types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->toc() << abort(FatalError);
    }

    return autoPtr<pdf>(cstrIter()(dict, rndGen));

    //return autoPtr<pdf>(new pdf);
}

// ************************************************************************* //

}
