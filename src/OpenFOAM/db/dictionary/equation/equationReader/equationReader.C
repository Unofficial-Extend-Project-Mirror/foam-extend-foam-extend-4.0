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
#include "mathematicalConstants.H"
#include "dimensionedScalar.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//const char* const Foam::equationReader::typeName = "equationReader";
    defineTypeNameAndDebug(equationReader, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::equationReader::fatalParseError
(
    const label index,
    const tokenList& tl,
    const label fromToken,
    const label toToken,
    const string& errorIn,
    const OStringStream& description
)
{
    OStringStream errorMessage;
    forAll(tl, i)
    {
        if (i == fromToken)
        {
            errorMessage << "<";
        }
        if (tl[i].isPunctuation() && tl[i].pToken() == token::COLON)
        {
            errorMessage << "^";
        }
        else
        {
            errorMessage << tl[i];
        }
        if (i == toToken)
        {
            errorMessage << ">";
        }
    }

    FatalErrorIn(errorIn) << "Parsing error in the equation for "
        << eqns_[index].equationName() << ", given by:" << endl
        << endl << token::TAB << eqns_[index].rawText() << endl
        << endl << "Error occurs withing the < angle brackets >:" << endl
        << endl << token::TAB << errorMessage.str() << endl << endl
        << description.str()
        << abort(FatalError);
}


Foam::string Foam::equationReader::stringPreconditioner(const string& rawText)
{
    string rawWorking(rawText);

    // Negative exponent workaround
    for (label i = 0; i < 10; i++)
    {
        string strTemp(name(i));
        string strTemp2(name(i));
        strTemp.append("e-");
        strTemp2.append("&");
        stringReplaceAll(rawWorking, strTemp, strTemp2);
    }
    
    stringReplaceAll(rawWorking, "^", " : ");
    stringReplaceAll(rawWorking, "(", " ( ");
    stringReplaceAll(rawWorking, ")", " ) ");
    stringReplaceAll(rawWorking, "+", " + ");
    stringReplaceAll(rawWorking, "-", " - ");
    stringReplaceAll(rawWorking, "&", "e-");
    stringReplaceAll(rawWorking, "*", " * ");
    stringReplaceAll(rawWorking, "/", " / ");
    stringReplaceAll(rawWorking, ",", " , ");
    
    // Leading negative workaround
    // This solution leads to dimensional problems
/*    IStringStream rawStream(rawWorking);
    token firstToken(rawStream);
    if (firstToken.isPunctuation() && firstToken.pToken() == token::SUBTRACT)
    {
        rawWorking.insert(0, "0");
    }
*/
    return rawWorking;
}


void Foam::equationReader::stringReplaceAll
(
    string& working,
    const string& findMe,
    const string& replaceWith
)
{
    size_t position(0);
    size_t offset(0);
    string subWorking(working);
    while (position != string::npos)
    {
        position = subWorking.string::find(findMe);
        if (position != string::npos)
        {
            working.::std::string::replace
            (
                position + offset, findMe.size(), replaceWith
            );
            offset += position + replaceWith.size();
            subWorking = subWorking.substr(position + findMe.size());
        }
    }
}


Foam::labelList Foam::equationReader::findMaxParenthesis
(
    const labelList& parenthesisList,
    const labelList& equationIndices
) const
{
    labelList returnMe(equationIndices.size());
    label currentMax(-1);
    label atIndex(0);
    bool groupDone(false);

    forAll(equationIndices, i)
    {
        if (mag(parenthesisList[equationIndices[i]]) > currentMax)
        {
            groupDone = false;
            atIndex = 0;
            currentMax = mag(parenthesisList[equationIndices[i]]);
            returnMe[atIndex] = equationIndices[i];
        }
        else if
        (
            (mag(parenthesisList[equationIndices[i]]) == currentMax)
         && (!groupDone)
        )
        {
            atIndex++;
            returnMe[atIndex] = equationIndices[i];
        }
        else if (mag(parenthesisList[equationIndices[i]]) < currentMax)
        {
            groupDone = true;
        }
    }
    returnMe.setSize(atIndex + 1);
    return returnMe;
}


Foam::labelList Foam::equationReader::findMaxOperation
(
    const labelList& opLvl,
    const labelList& indices
)
{
    label maxOpLvl(-1);
    label atIndex(-1);
    labelList returnMe(indices.size());
    bool groupDone(false);

    forAll(indices, i)
    {
        if (opLvl[indices[i]] > maxOpLvl)
        {
            atIndex = 0;
            returnMe[0] = indices[i];
            if (i > 0)
            {
                atIndex++;
                returnMe[0] = indices[i - 1];
                returnMe[1] = indices[i];
            }
            groupDone = false;
            maxOpLvl = opLvl[indices[i]];
        }
        else if
        (
            (
                (opLvl[indices[i]] == maxOpLvl)
             || (opLvl[indices[i]] == 0)
            )
         && !groupDone
        )
        {
            atIndex++;
            returnMe[atIndex] = indices[i];
        }
        else if ((opLvl[indices[i]] < maxOpLvl) && (opLvl[indices[i]] != 0))
        {
            groupDone = true;
        }
    }
    returnMe.setSize(atIndex + 1);
    return returnMe;
}


void Foam::equationReader::absorbNegatives
(
    const label equationIndex,
    const tokenList& tl,
    labelList& eqnIndices,
    labelList& subEqnIndices,
    equationOperationList& map,
    labelList& opLvl
)
{
    // Negatives are identified by a map with a negative sourceIndex
    // To accommodate this behaviour, source indices are 1-indexed in the map.
    forAll(subEqnIndices, i)
    {
        if (map[subEqnIndices[i]].dictLookupIndex() == -1)
        {
            if
            (
                (subEqnIndices.size() == i)
             || (opLvl[subEqnIndices[i + 1]] != 0)
            )
            {
                OStringStream description;
                description << "Misplaced negative / subtraction operator.";
                fatalParseError
                (
                    equationIndex,
                    tl,
                    subEqnIndices[i],
                    subEqnIndices[i],
                    "equationReader::absorbNegatives",
                    description
                );
            }
            map[subEqnIndices[i + 1]].sourceIndex() = 
                -map[subEqnIndices[i + 1]].sourceIndex();

            trimListWithParent(eqnIndices, subEqnIndices, i, i);
        }
    }
}


void Foam::equationReader::dsEqual
(
    dimensionedScalar& dso,
    const dimensionedScalar& dsi
)
{
    dso.value() = dsi.value();
    dso.name() = dsi.name();
    dso.dimensions().reset(dsi.dimensions());
}


void Foam::equationReader::trimListWithParent
(
    labelList& parent,
    labelList& indices,
    label from,
    label to,
    label exceptFor
)
{
    if
    (
        (to > (indices.size() - 1))
     || (from < 0)
     || (from > (indices.size() - 1))
     || (to < from)
    )
    {
        FatalErrorIn
        (
            "equationReader::trimListWithParent(parent, indices, from, to, "
            "exceptFor)"
        )
            << "Bad indices.  parent is " << parent << ", indices are "
            << indices << ", from is " << from << ", to is " << to
            << " exceptFor is " << exceptFor << "."
            << abort(FatalError);
    }

    for (label i(from); i <= to; i++)
    {
        label removeMe(indices[i]);
        if (i == exceptFor) continue;
        forAll(parent, j)
        {
            if (parent[j] == removeMe)
            {
                trimList(parent, j, j);
                break;
            }
        }
    }
    trimList(indices, from, to, exceptFor);
}


void Foam::equationReader::trimList
(
    labelList& indices,
    label from,
    label to,
    label exceptFor
)
{
    if
    (
        (to > (indices.size() - 1))
     || (from < 0)
     || (from > (indices.size() - 1))
     || (to < from)
    )
    {
        FatalErrorIn
        (
            "equationReader::trimList(indices, from, to, exceptFor"
            "exceptFor)"
        )
            << "Bad indices.  indices are " << indices << ", from is "
            << from << ", to is " << to << " exceptFor is " << exceptFor << "."
            << abort(FatalError);
    }

    if (!(exceptFor == from && from == to))
    {
        
        if (exceptFor == from)
        {
            from++;
        }
        else if (exceptFor == to)
        {
            to--;
        }
        else if ((exceptFor < to) && (exceptFor > from))
        {
            indices[from] = indices[exceptFor];
            from++;
        }

        for (label i(from); i < (indices.size() + from - to - 1); i++)
        {
            indices[i] = indices[i + to - from + 1];
        }
        indices.setSize(indices.size() - to + from - 1);
    }
}


Foam::label Foam::equationReader::findIndex
(
    const label value,
    const labelList& indexList
) const
{
    forAll (indexList, i)
    {
        if (indexList[i] == value)
        {
            return i;
        }
    }
    return -1;
}


Foam::equationOperation Foam::equationReader::findSource
(
    const word& varName
)
{
    // Search order:
    // -other equations
    // -externalDScalars
    // -externalScalars
    // -externalScalarLists
    // -dictSources

    // Searching known equations
    for (label eqs(0); eqs < eqns_.size(); eqs++)
    {
        if (eqns_[eqs].equationName() == varName)
        {
            return equationOperation
            (
                equationOperation::slequation,
                eqs + 1,
                0,
                equationOperation::otnone
            );
        }
    }
    
    forAll(externalDScalars_, i)
    {
        if (externalDScalars_[i].name() == varName)
        {
            return equationOperation
            (
                equationOperation::slexternalDScalar,
                i + 1,
                0,
                equationOperation::otnone
            );
        }
    }
    
    if (externalScalars_.size() != externalScalarNames_.size())
    {
        FatalErrorIn("equationReader::findSource")
            << "Size mismatch detected in externalScalars = "
            << externalScalars_.size() << " and externalScalarNames = "
            << externalScalarNames_.size() << "."
            << abort(FatalError);
    }
    
    forAll(externalScalars_, j)
    {
        if (externalScalarNames_[j] == varName)
        {
            return equationOperation
            (
                equationOperation::slexternalScalar,
                j + 1,
                0,
                equationOperation::otnone
            );
        }
    }
    
    forAll(externalScalarListNames_, i)
    {
        if (externalScalarListNames_[i] == varName)
        {
            return equationOperation
            (
                equationOperation::slexternalScalarList,
                i + 1,
                0,
                equationOperation::otnone
            );
        }
    }

    forAll(dictSources_, i)
    {
        if (dictSources_[i].found(varName) && !dictSources_[i].isDict(varName))
        {
            label dictLookupIndex(-1);
//            wordList& wl(eqns_[equationIndex].dictLookups());

            ITstream srcStrm
            (
                dictSources_[i].lookup(varName)
            );
            if (isDimensionedScalar(srcStrm) || isScalar(srcStrm))
            {
                // Look for varName in dictLoookup names, append if not found
                forAll(dictLookups_, j)
                {
                    if (varName == dictLookups_[j])
                    {
                        dictLookupIndex = j;
                        break;
                    }
                }
                if (dictLookupIndex < 0)
                {
                    dictLookups_.setSize(dictLookups_.size() + 1);
                    dictLookups_.set
                    (
                        dictLookups_.size() - 1,
                        new word(varName)
                    );
                    dictLookupIndex = dictLookups_.size() - 1;
                }
                
                return equationOperation
                (
                    equationOperation::sldictSource,
                    i + 1,
                    dictLookupIndex,
                    equationOperation::otnone
                );
            }
            else if (isEquation(srcStrm))
            {
                // Is it a known equation already?
                for (label eqs(0); eqs < eqns_.size(); eqs++)
                {
                    if (eqns_[eqs].equationName() == varName)
                    {
                        return equationOperation
                        (
                            equationOperation::slequation,
                            eqs + 1,
                            0,
                            equationOperation::otnone
                        );
                    }
                }
                // The equation has not been read yet.  Create an unparsed
                // equation.  It will be parsed during evaluate, or if it
                // is later read.
                equation eqn(srcStrm);
                eqn.equationName() = varName;
                createEquation(eqn);
/*
                token it(srcStrm);
                createEquation
                (
                    equation
                    (
                        varName,
                        it.stringToken()
                    )
                );
*/
                return equationOperation
                (
                    equationOperation::slequation,
                    eqns_.size(),
                    0,
                    equationOperation::otnone
                );
            }
        } // end if varName is in dictionary
    } // end dictionary search loop
    return equationOperation
    (
        equationOperation::slnone,
        0,
        0,
        equationOperation::otnone
    );

}


Foam::label Foam::equationReader::addInternalScalar(const scalar& value)
{
    forAll(internalScalars_, i)
    {
        if (mag(internalScalars_[i] - value) < VSMALL)
        {
            return i;
        }
    }
    internalScalars_.setSize(internalScalars_.size() + 1);
    internalScalars_.set(internalScalars_.size() - 1, new scalar(value));
    return internalScalars_.size() - 1;
}


bool Foam::equationReader::isDimensionedScalar(ITstream& is)
{
    tokenList tl(12);
    label found(0);
    while (!is.eof())
    {
        tl[found] = token(is);
        found++;
        if (found > 12)
        {
            is.rewind();
            return false;
        }
    }
    if
    (
        (
            (found == 11)
         && tl[0].isWord()
         && tl[1].isPunctuation()
         && tl[2].isNumber()
         && tl[3].isNumber()
         && tl[4].isNumber()
         && tl[5].isNumber()
         && tl[6].isNumber()
         && tl[7].isNumber()
         && tl[8].isNumber()
         && tl[9].isPunctuation()
         && tl[10].isNumber()
        )
     || (
            (found == 9)
         && tl[0].isWord()
         && tl[1].isPunctuation()
         && tl[2].isNumber()
         && tl[3].isNumber()
         && tl[4].isNumber()
         && tl[5].isNumber()
         && tl[6].isNumber()
         && tl[7].isPunctuation()
         && tl[8].isNumber()
        )
    )
    {
        is.rewind();
        return true;
    }
    else
    {
        is.rewind();
        return false;
    }
}


bool Foam::equationReader::isScalar(ITstream& is)
{
    token firstToken(is);
    if (firstToken.isNumber() && is.eof())
    {
        is.putBack(firstToken);
        return true;
    }
    else
    {
        is.putBack(firstToken);
        return false;
    }
}


bool Foam::equationReader::isEquation(ITstream& is)
{
    token firstToken(is);
    if (firstToken.isString() && is.eof())
    {
        is.putBack(firstToken);
        return true;
    }
    else
    {
        is.putBack(firstToken);
        tokenList tl(12);
        label found(0);
        while (!is.eof())
        {
            tl[found] = token(is);
            found++;
            if (found == 12)
            {
                is.rewind();
                return false;
            }
        }
        if
        (
            (
                (found == 11)
             && tl[0].isWord()
             && tl[1].isPunctuation()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isNumber()
             && tl[7].isNumber()
             && tl[8].isNumber()
             && tl[9].isPunctuation()
             && (tl[10].isString() || tl[10].isNumber())
            )
         || (
                (found == 10)
             && tl[0].isPunctuation()
             && tl[1].isNumber()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isNumber()
             && tl[7].isNumber()
             && tl[8].isPunctuation()
             && (tl[9].isString() || tl[9].isNumber())
            )
         || (
                (found == 9)
             && tl[0].isWord()
             && tl[1].isPunctuation()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isNumber()
             && tl[7].isPunctuation()
             && (tl[8].isString() || tl[8].isNumber())
            )
         || (
                (found == 8)
             && tl[0].isPunctuation()
             && tl[1].isNumber()
             && tl[2].isNumber()
             && tl[3].isNumber()
             && tl[4].isNumber()
             && tl[5].isNumber()
             && tl[6].isPunctuation()
             && (tl[7].isString() || tl[8].isNumber())
            )
        )
        {
            is.rewind();
            return true;
        }
        else
        {
            is.rewind();
            return false;
        }
    }
}


void Foam::equationReader::reportOperationDisabled
(
    const label& index,
    const label& i,
    const dimensionedScalar& ds
) const
{
    // do nothing
}


void Foam::equationReader::reportOperationEnabled
(
    const label& index,
    const label& i,
    const dimensionedScalar& ds
) const
{
    Info << "Performing operation: ["
        << equationOperation::opName
        (
            eqns_[index].ops()[i].operation()
        ) << "] using source [";
    if
    (
        eqns_[index].ops()[i].sourceList()
     == equationOperation::sldictSource
    )
    {
        Info << dictLookups_[eqns_[index].ops()[i].dictLookupIndex()];
    }
    else if
    (
        eqns_[index].ops()[i].sourceList()
     == equationOperation::slequation
    )
    {
        Info << eqns_
            [
                mag(eqns_[index].ops()[i].sourceIndex() - 1)
            ].equationName();
    }
    else if
    (
        eqns_[index].ops()[i].sourceList()
     == equationOperation::slexternalScalar
    )
    {
        Info << externalScalarNames_
            [
                mag(eqns_[index].ops()[i].sourceIndex() - 1)
            ];
    }
    else
    {
        Info << ds.name();
    }
    Info << "] read from ["
        << equationOperation::sourceName
        (
            eqns_[index].ops()[i].sourceList()
        ) << "]..." << endl;
}


void Foam::equationReader::reportResultDisabled
(
    const dimensionedScalar& ds
) const
{
    // do nothing
}


void Foam::equationReader::reportResultEnabled
(
    const dimensionedScalar& ds
) const
{
    Info << "Operaion result is " << ds << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationReader::equationReader()
{
    if (debug > 1)
    {
        reportOperationFunction_
            = &Foam::equationReader::reportOperationEnabled;
        reportResultFunction_
            = &Foam::equationReader::reportResultEnabled;
    }
    else
    {
        reportOperationFunction_
            = &Foam::equationReader::reportOperationDisabled;
        reportResultFunction_
            = &Foam::equationReader::reportResultDisabled;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationReader::~equationReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::equationReader::addDataSource(const dictionary& dict)
{
    dictSources_.setSize(dictSources_.size() + 1);
    dictSources_.set(dictSources_.size() - 1, &dict);
}


void Foam::equationReader::addDataSource
(
    const scalar& value,
    const word& name,
    const dimensionSet& dimensions
)
{
    equationOperation source(findSource(name));
    if
    (
        source.sourceList() != equationOperation::slnone
     && source.sourceList() != equationOperation::sldictSource
    )
    {
        FatalErrorIn("equationReader::addDataSource")
            << "Attempting to add external scalar named " << name << " when "
            << "this name already exists as a "
            << equationOperation::sourceName(source.sourceList()) << " source."
            << abort(FatalError);
    }
    
    label newSize(externalScalars_.size() + 1);
    
    externalScalars_.setSize(newSize);
    externalScalarNames_.setSize(newSize);
    externalScalarDimensions_.setSize(newSize);
    
    externalScalars_.set(newSize - 1, &value);
    externalScalarNames_[newSize - 1] = name;
    externalScalarDimensions_.set(newSize - 1, new dimensionSet(dimensions));
}


void Foam::equationReader::addDataSource(const dimensionedScalar& ds)
{
    addDataSource(ds.value(), ds.name(), ds.dimensions());
}


void Foam::equationReader::addDataSource
(
    const scalarList& slist,
    const word& name,
    const dimensionSet& dimensions
)
{
    equationOperation source(findSource(name));
    if
    (
        source.sourceList() != equationOperation::slnone
     && source.sourceList() != equationOperation::sldictSource
    )
    {
        FatalErrorIn("equationReader::addDataSource")
            << "Attempting to add external scalarList named " << name
            << " when this name already exists as a "
            << equationOperation::sourceName(source.sourceList()) << " source."
            << abort(FatalError);
    }

    label newSize(externalScalarLists_.size() + 1);
    
    externalScalarLists_.setSize(newSize);
    externalScalarListNames_.setSize(newSize);
    externalScalarListDimensions_.setSize(newSize);
    externalScalarListIndex_.setSize(newSize);

    externalScalarLists_.set(newSize - 1, &slist);
    externalScalarListNames_[newSize - 1] = name;
    externalScalarListDimensions_.set
    (
        newSize - 1,
        new dimensionSet(dimensions)
    );
    externalScalarListIndex_[newSize - 1] = 0;    
}


void Foam::equationReader::setListIndex(const word& name, label newIndex)
{
    equationOperation source(findSource(name));
    if (source.sourceList() != equationOperation::slexternalScalarList)
    {
        FatalErrorIn("equationReader::setListIndex")
            << "setListIndex called for variable " << name << " which is not "
            << "a scalarList source.  Its source is: "
            << equationOperation::sourceName(source.sourceList())
            << abort(FatalError);
    }
    label maxIndex(externalScalarLists_[source.sourceIndex() - 1].size());
    if (newIndex < 0 || newIndex > maxIndex)
    {
        FatalErrorIn("equationReader::setListIndex")
            << "Index " << newIndex << " is out of bounds (0, " << maxIndex - 1
            << ") for scalarList " << name << "."
            << abort(FatalError);
    }
    externalScalarListIndex_[source.sourceIndex() - 1] = newIndex;
}


void Foam::equationReader::setListIndex(label newIndex)
{
    forAll(externalScalarListNames_, i)
    {
        setListIndex(externalScalarListNames_[i], newIndex);
    }
}


void Foam::equationReader::createEquation
(
    equation eqn
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        if (debug)
        {
            Info << "Creating equation " << eqn.equationName() << " at index "
                << eqns_.size() << endl;
        }
        label newSize(eqns_.size() + 1);
        
        eqns_.setSize(newSize);
        outputScalars_.setSize(newSize);
        outputScalarNames_.setSize(newSize);
        outputScalarDimensions_.setSize(newSize);
        outputScalarLists_.setSize(newSize);
        outputScalarListDimensions_.setSize(newSize);

        eqns_.set(newSize - 1, new equation(eqn));
        outputScalars_.set(newSize - 1, NULL);
        outputScalarNames_[newSize - 1] = word::null;
        outputScalarDimensions_.set
        (
            newSize - 1,
            new dimensionSet(dimless)
        );
        outputScalarLists_.set(newSize - 1, NULL);
        outputScalarListDimensions_.set
        (
            newSize - 1,
            new dimensionSet(dimless)
        );
    }
    else
    {
        FatalErrorIn("equationReader::createEquation")
            << "Equation " << eqn.equationName() << " already exists."
            << abort(FatalError);
    }
}


void Foam::equationReader::linkOutput
(
    const word& eqnName,
    dimensionedScalar& outputDVar
)
{
    label index(lookup(eqnName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::linkOutput")
            << "Equation name " << eqnName << "not found."
            << abort(FatalError);
    }
    linkOutput
    (
        index,
        outputDVar.value(),
        outputDVar.name(),
        outputDVar.dimensions()
    );
}


void Foam::equationReader::linkOutput
(
    const word& eqnName,
    scalar& value,
    const word& name,
    const dimensionSet& dimensions
)
{
    label index(lookup(eqnName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::linkOutput")
            << "Equation name " << eqnName << "not found."
            << abort(FatalError);
    }
    linkOutput
    (
        index,
        value,
        name,
        dimensions
    );
}


void Foam::equationReader::linkOutput
(
    const word& eqnName,
    scalarList& outputSList,
    const dimensionSet& dimensions
)
{
    label index(lookup(eqnName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::linkOutput")
            << "Equation name " << eqnName << "not found."
            << abort(FatalError);
    }
    linkOutput(index, outputSList, dimensions);
}


void Foam::equationReader::linkOutput
(
    label index,
    dimensionedScalar& outputDVar
)
{
    linkOutput
    (
        index,
        outputDVar.value(),
        outputDVar.name(),
        outputDVar.dimensions()
    );
}


void Foam::equationReader::linkOutput
(
    label index,
    scalar& value,
    const word& name,
    const dimensionSet& dimensions
)
{
    if (index < 0 || index >= eqns_.size())
    {
        FatalErrorIn("equationReader::linkOutput")
            << "Index " << index << " out of range (0, " << eqns_.size() - 1
            << ")."
            << abort(FatalError);
    }
    outputScalars_.set(index, &value);
    outputScalarNames_[index] = name;
    outputScalarDimensions_[index].reset(dimensions);
    outputScalarLists_.set(index, NULL);
}


void Foam::equationReader::linkOutput
(
    label index,
    scalarList& outputSList,
    const dimensionSet& dimensions
)
{
    outputScalarLists_.set(index, & outputSList);
    outputScalarListDimensions_.set(index, new dimensionSet(dimensions));
    outputScalars_.set(index, NULL);
}


void Foam::equationReader::readEquation
(
    equation eqn
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        createEquation(eqn);
    }
    else if(eqn.rawText() != eqns_[index].rawText())
    {
        eqns_[index].clear();
    }
}


void Foam::equationReader::readEquation
(
    equation eqn,
    dimensionedScalar& outputDVar
)
{
    readEquation
    (
        eqn,
        outputDVar.value(),
        outputDVar.name(),
        outputDVar.dimensions()
    );
}


void Foam::equationReader::readEquation
(
    equation eqn,
    scalar& value,
    const word& name,
    const dimensionSet& dimensions
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        index = eqns_.size();
        createEquation(eqn);
    }
    else if(eqn.rawText() != eqns_[index].rawText())
    {
        eqns_[index].clear();
    }
    linkOutput(index, value, name, dimensions);
}


void Foam::equationReader::readEquation
(
    equation eqn,
    scalarList& outputSList,
    const dimensionSet& dimensions
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        index = eqns_.size();
        createEquation(eqn);
    }
    else if(eqn.rawText() != eqns_[index].rawText())
    {
        eqns_[index].clear();
    }
    linkOutput(index, outputSList, dimensions);
}


void Foam::equationReader::readEquation
(
    const dictionary& sourceDict,
    const word& eqnName
)
{
    equation eqn(sourceDict.lookup(eqnName));
    eqn.equationName() = eqnName;
    readEquation(eqn);
}


void Foam::equationReader::readEquation
(
    const dictionary& sourceDict,
    const word& eqnName,
    dimensionedScalar& outputDVar
)
{
    equation eqn(sourceDict.lookup(eqnName));
    eqn.equationName() = eqnName;
    readEquation
    (
        eqn,
        outputDVar.value(),
        outputDVar.name(),
        outputDVar.dimensions()
    );
}


void Foam::equationReader::readEquation
(
    const dictionary& sourceDict,
    const word& eqnName,
    scalar& value,
    const word& name,
    const dimensionSet& dimensions
)
{
    equation eqn(sourceDict.lookup(eqnName));
    eqn.equationName() = eqnName;
    readEquation
    (
        eqn,
        value,
        name,
        dimensions
    );
}

/*
void Foam::equationReader::createEquation
(
    equation eqn,
    dimensionedScalar * outputDVar
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        if (debug)
        {
            Info << "Creating equation " << eqn.equationName() << " at index "
                << eqns_.size() << endl;
        }
        eqns_.setSize(eqns_.size() + 1);
        eqns_.set(eqns_.size() - 1, new equation(eqn));

        outputScalars_.setSize(eqns_.size());
        outputScalarNames_.setSize(eqns_.size());
        outputScalarDimensions_.setSize(eqns_.size());
        if (outputDVar != NULL)
        {
            outputScalars_.set(eqns_.size() - 1, & outputDVar->value());
            outputScalarNames_[eqns_.size() - 1] = outputDVar->name();
            outputScalarDimensions_.set
            (
                eqns_.size() - 1,
                new dimensionSet(outputDVar->dimensions())
            );
        }
        else
        {
            outputScalars_.set(eqns_.size() - 1, NULL);
            outputScalarNames_[eqns_.size() - 1] = word::null;
            outputScalarDimensions_.set
            (
                eqns_.size() - 1,
                new dimensionSet(dimless)
            );
        }
    }
    else
    {
        FatalErrorIn("equationReader::createEquation")
            << "Equation " << eqn.equationName() << " already exists.  Use "
            << "readEquation to suppress error."
            << abort(FatalError);
    }
}

void Foam::equationReader::createEquation
(
    equation eqn,
    scalar * outputVar,
    const word& outputName,
    const dimensionSet& outputDimensions
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        if (debug)
        {
            Info << "Creating equation " << eqn.equationName() << " at index "
                << eqns_.size() << endl;
        }
        eqns_.setSize(eqns_.size() + 1);
        eqns_.set(eqns_.size() - 1, new equation(eqn));

        outputScalars_.setSize(eqns_.size());
        outputScalarNames_.setSize(eqns_.size());
        outputScalarDimensions_.setSize(eqns_.size());

        outputScalars_.set(eqns_.size() - 1, outputVar);
        outputScalarNames_[eqns_.size() - 1] = outputName;
        outputScalarDimensions_.set
        (
            eqns_.size() - 1,
            new dimensionSet(outputDimensions)
        );
    }
    else
    {
        FatalErrorIn("equationReader::createEquation")
            << "Equation already exists.  Use readEquation to suppress error."
            << abort(FatalError);
    }
}


void Foam::equationReader::readEquation
(
    const dictionary& sourceDict,
    const word& eqnName,
    dimensionedScalar * outputDVar
)
{
    equation eqn(sourceDict.lookup(eqnName));
    eqn.equationName() = eqnName;

    readEquation
    (
        eqn,
        outputDVar
    );
}


void Foam::equationReader::readEquation
(
    const dictionary& sourceDict,
    const word& eqnName,
    scalar * outputVar,
    const word& outputName,
    const dimensionSet& outputDimensions
)
{
    equation eqn(sourceDict.lookup(eqnName));
    eqn.equationName() = eqnName;

    readEquation
    (
        eqn,
        outputVar,
        outputName,
        outputDimensions
    );
}


void Foam::equationReader::readEquation
(
    equation eqn,
    dimensionedScalar * outputDVar
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        createEquation(eqn, outputDVar);
        parse(eqns_.size() - 1);
    }
    else if
    (
        (eqn.rawText() != eqns_[index].rawText())
     || (!eqns_[index].size())
    )
    {
        eqns_[index].rawText() = eqn.rawText();
        parse(index);
        if (outputDVar != NULL)
        {
            outputScalars_.set(index, & outputDVar->value());
            outputScalarNames_[index] = outputDVar->name();
            outputScalarDimensions_.set
            (
                index,
                new dimensionSet(outputDVar->dimensions())
            );
        }
    }
}


void Foam::equationReader::readEquation
(
    equation eqn,
    scalar * outputVar,
    const word& outputName,
    const dimensionSet& outputDimensions
)
{
    label index(lookup(eqn.equationName()));
    if (index < 0)
    {
        createEquation(eqn, outputVar, outputName, outputDimensions);
        parse(eqns_.size() - 1);
    }
    else if
    (
        (eqn.rawText() != eqns_[index].rawText())
     || (!eqns_[index].size())
    )
    {
        parse(index);
        outputScalars_.set(index, outputVar);
        outputScalarNames_[index] = outputName;
        outputScalarDimensions_.set
        (
            index,
            new dimensionSet(outputDimensions)
        );
    }
}
*/


void Foam::equationReader::update(const word& equationName)
{
    label index(lookup(equationName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::update(const word)")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    update(index);
}


void Foam::equationReader::update(const label& index)
{
    if (outputScalars_.set(index))
    {
        dimensionedScalar eval = evaluate(index);
        outputScalars_[index] = eval.value();
        
        if
        (
            !eqns_[index].changeDimensions()
         && dimensionSet::debug
         && outputScalarDimensions_[index] != eval.dimensions()
        )
        {
            WarningIn("equationReader::update")
                << "Dimension error thrown for equation "
                << eqns_[index].equationName() << ", given by:"
                << token::NL << token::TAB
                << eqns_[index].rawText();

            outputScalarDimensions_[index] = eval.dimensions();
        }
    }
    else if (outputScalarLists_.set(index))
    {
        evaluateField
        (
            index,
            outputScalarLists_[index],
            outputScalarDimensions_[index]
        );
    }
}


void Foam::equationReader::update()
{
    for (label i = 0; i < eqns_.size(); i++)
    {
        update(i);
    }
}


bool Foam::equationReader::found(const word& equationName)
{
    for (label i = 0; i < eqns_.size(); i++)
    {
        if (eqns_[i].equationName() == equationName)
        {
            return true;
        }
    }
    return false;
}


Foam::label Foam::equationReader::lookup(const word& equationName)
{
    for (label i = 0; i < eqns_.size(); i++)
    {
        if (eqns_[i].equationName() == equationName)
        {
            return i;
        }
    }

    return -1;
}


void Foam::equationReader::deleteEquation(const word& equationName)
{
    label index(lookup(equationName));
    if (index < 0)
    {
        WarningIn("equationReader::deleteEquation(equationName)")
            << "Equation name " << equationName << " not found." << endl;
    }
    deleteEquation(index);
}


void Foam::equationReader::deleteEquation(const label& index)
{
    if ((index < 0) || (index >= eqns_.size()))
    {
        FatalErrorIn("equationReader::deleteEquation(index)")
            << "Index " << index << " out of bounds (0, " << eqns_.size() - 1
            << ")"
            << abort(FatalError);
    }
    for (label i = index; i < (eqns_.size() - 1); i++)
    {
        eqns_[i] = eqns_[i + 1];
    }
    
    eqns_.setSize(eqns_.size() - 1);
}


/*
void Foam::equationReader::status()
{
    Info << "*** equationReader object status report ***" << endl;
    Info << "Data sources:" << endl;
    Info << token::TAB << dictSources_.size() << " dictionaries, with recorded "
        << "keywords:" << endl;
    forAll(dictLookups_, i)
    {
        Info << token::TAB << token::TAB << dictLookups_[i] << endl;
    }
    Info << token::TAB << externalDScalars_.size() << " external "
        << "dimensionedScalars with names:" << endl;
    forAll(externalDScalars_, i)
    {
        Info << token::TAB << token::TAB << externalDScalars_[i].name()
            << endl;
    }
    Info << token::TAB << externalScalars_.size() << " external "
        << "scalars with names:" << endl;
    forAll(externalScalarNames_, i)
    {
        Info << token::TAB << token::TAB << externalScalarNames_[i] << endl;
    }
    Info << token::TAB << internalScalars_.size() << " internal scalars."
        << endl;
    Info << eqns_.size() << " equations in memory:" << endl;
    forAll(eqns_, i)
    {
        Info << i << ":" << token::TAB << eqns_[i].size() << " operations,"
            << token::TAB << eqns_[i].equationName() << token::TAB;
        if ((outputScalars_.size() > i) && outputScalars_.set(i))
        {
            Info << "active output to " << outputScalarNames_[i] << "."
                << endl;
        }
        else
        {
            Info << "passive output." << endl;
        }
        if (eqns_[i].changeDimensions())
        {
            Info << token::TAB << "Dimension override to: "
                << eqns_[i].overrideDimensions() << endl;
        }
        Info << token::TAB << eqns_[i].rawText() << endl;
    }
}
*/

/*
template <>
dimensioned<scalar>::dimensioned
(
    Istream& is
)
:
    name_(is),
    dimensions_(dimless),
    value_(0)
{
    is.rewind();
    token t(is);

    if (t.isString())
    {
        equationReader eqn;
        eqn.readEquation
        (
            equation
            (
                name_,
                t.stringToken()
            )
        );
        *this = eqn.evaluate(0);
    }
    else
    {
        dimensions_ = dimensionSet(is);
        value_ = scalar(pTraits<scalar>(is));
    }
}
*/

#include "equationReaderAssignPointerFunctions.C"
#include "equationReaderCreateMap.C"
#include "equationReaderEvaluate.C"
#include "equationReaderGetSource.C"
#include "equationReaderParse.C"

// ************************************************************************* //
