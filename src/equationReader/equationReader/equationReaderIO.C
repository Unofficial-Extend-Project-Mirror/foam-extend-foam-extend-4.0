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

#include "equation.H"
#include "equationOperation.H"
#include "equationReader.H"
#include "fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::equationReader::dataSourceStatus(Ostream& os) const
{
    dictionary dict;
    dict.set("activeSources", activeSourceNames_);
    fileNameList dictPaths(dictSources_.size());
    forAll(dictSources_, i)
    {
        dictPaths[i] = dictSources_[i].name();
    }
    dict.set("dictSources", dictPaths);
    dict.set("dictLookups", dictLookups_);
    dict.set("scalarSources", scalarSources_.outputDictionary());
    dict.set("vectorSources", vectorSources_.outputDictionary());
    dict.set("tensorSources", tensorSources_.outputDictionary());
    dict.set("diagTensorSources", diagTensorSources_.outputDictionary());
    dict.set("symmTensorSources", symmTensorSources_.outputDictionary());
    dict.set
    (
        "sphericalTensorSources",
        sphericalTensorSources_.outputDictionary()
    );
    dict.set("cellIndex", cellIndex_);
    dict.set("internalScalars", internalScalars_);
    dictionary superDict;
    superDict.set("dataSources", dict);
    os << superDict;
    return os;
}


// * * * * * * * * * * * * * Friend IOstream Operators * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, equationReader& I)
{
    dictionary dict(is);
    dictionary eqnsDict(dict.subDict("equations"));
    wordList eqnNames(eqnsDict.toc());
    forAll(eqnNames, i)
    {
        equation eq
        (
            eqnsDict.subDict(eqnNames[i]).lookup(eqnNames[i]),
            eqnNames[i]
        );
        I.createEquation(eq);
    }
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const equationReader& I)
{
    dictionary eqnsDict;
    forAll(I, i)
    {
        dictionary eqnEntry;
        eqnEntry.set(I[i].name(), I[i]);
        eqnEntry.set("lastResult", I[i].lastResult());
        eqnsDict.set(I[i].name(), eqnEntry);
    }
    dictionary dict;
    dict.set("equations", eqnsDict);
    os << dict;
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
