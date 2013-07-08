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

#include "equationReader.H"
#include "fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::equationReader::dataSourceStatus(Ostream& os) const
{
    dictionary dict;
    fileNameList dictPaths(dictSources_.size());
    forAll(dictSources_, i)
    {
        dictPaths[i] = dictSources_[i].name();
    }
    dict.set("dictSources", dictPaths);
    dict.set("dictLookups", dictLookups_);
    dict.set("externalDScalars", externalDScalars_);
    dict.set("externalScalars", externalScalars_);
    dict.set("externalScalarNames", externalScalarNames_);
    dict.set("externalScalarListName", externalScalarListNames_);
    dict.set("externalScalarListDimensions", externalScalarListDimensions_);
    dict.set("externalScalarListIndices", externalScalarListIndex_);
    PtrList<dimensionedScalar> outputScalarList;
    forAll(eqns_, i)
    {
        if (outputScalars_.set(i))
        {
            outputScalarList.setSize(outputScalarList.size() + 1);
            outputScalarList.set
            (
                outputScalarList.size() - 1,
                new dimensionedScalar
                (
                    outputScalarNames_[i],
                    outputScalarDimensions_[i],
                    outputScalars_[i]
                )
            );
        }
    }
    dict.set("outputScalars", outputScalarList);
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
    forAll(I.eqns_, i)
    {
        dictionary eqnEntry;
        eqnEntry.set(I.eqns_[i].equationName(), I.eqns_[i]);
        eqnEntry.set("lastResult", I.eqns_[i].lastResult());
        if (I.outputScalars_.set(i))
        {
            dimensionedScalar ds
            (
                I.outputScalarNames_[i],
                I.outputScalarDimensions_[i],
                I.outputScalars_[i]
            );
            eqnEntry.set("outputVar", ds);
        }
        eqnsDict.set(I.eqns_[i].equationName(), eqnEntry);
    }
    dictionary dict;
    dict.set("equations", eqnsDict);
    os << dict;
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

