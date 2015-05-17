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

#include "objectRegistry.H"
#include "equationReader.H"
#include "mathematicalConstants.H"
#include "dimensionedScalar.H"
#include "equation.H"
#include "equationOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
//const char* const Foam::equationReader::typeName = "equationReader";
    defineTypeNameAndDebug(equationReader, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::equationReader::createEquation
(
    equation& eqn
) const
{
    label index(lookup(eqn.name()));
    if (index < 0)
    {
        if (debug)
        {
            Info << "Creating equation " << eqn.name() << " at index "
                << size() << endl;
        }
        label newSize(size() + 1);

        setSize(newSize);
        set(newSize - 1, new equation(eqn));

        evaluateDimsFunctions_.setSize(newSize);
        evaluateDimsFunctions_.set
        (
            newSize - 1,
            new evaluateDimsFunction
            (
                eqn.changeDimensions()
              ? &Foam::equationReader::evaluateDimsDisabled
              : &Foam::equationReader::evaluateDimsEnabled
            )
        );
        return newSize - 1;
    }
    else
    {
        FatalErrorIn("equationReader::createEquation")
            << "Equation " << eqn.name() << " already exists."
            << abort(FatalError);
    }
    return -1;
}


void Foam::equationReader::fatalParseError
(
    const label index,
    const tokenList& tl,
    const label fromToken,
    const label toToken,
    const string& errorIn,
    const OStringStream& description
) const
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
        << operator[](index).name() << ", given by:" << endl
        << endl << token::TAB << operator[](index).rawText() << endl
        << endl << "Error occurs withing the < angle brackets >:" << endl
        << endl << token::TAB << errorMessage.str() << endl << endl
        << description.str()
        << abort(FatalError);
}


Foam::string Foam::equationReader::stringPreconditioner
(
    const string& rawText
)
{
    string rawWorking(rawText);

    // Negative exponent workaround
    for (label i = 0; i < 10; i++)
    {
        string strTemp(name(i));
        string strTempU(name(i));
        string strTemp2(name(i));
        strTemp.append("e-");
        strTempU.append("E-");
        strTemp2.append("&");
        stringReplaceAll(rawWorking, strTemp, strTemp2);
        stringReplaceAll(rawWorking, strTempU, strTemp2);
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


void Foam::equationReader::removePowExponents
(
    const label index,
    tokenList& tl,
    PtrList<equationOperation>& map,
    labelList& opLvl,
    labelList& pl
) const
{
    // Remove pow(a,b) exponent part 'b' from an equation and create a sub-
    // equation.
    label tokenI(0);

    while (tokenI < map.size())
    {
        if (map[tokenI].operation() == equationOperation::otpow)
        {
            // Found a 'pow('. Look for ','; fail on ')', or end of list
            // pl checks ensure the ',' or ')' relate to the 'pow(', and not
            // another function / parethesis
            const label powFoundAt(tokenI);
            const label pLvl(pl[tokenI]);
            while ((opLvl[tokenI] != 5) || (pl[tokenI] != pLvl))
            {
                if
                (
                    ((opLvl[tokenI] == -4) && (pl[tokenI] == pLvl))
                 || (tokenI == (map.size() - 1))
                )
                {
                    OStringStream description;
                    description << "pow() function takes two arguments.";
                    fatalParseError
                    (
                        index,
                        tl,
                        powFoundAt,
                        tokenI,
                        "equationReader::removePowExponents",
                        description
                    );
                }
                tokenI++;
            }

            // Found 'pow( ... ,' look for ')', fail on list end
            const label commaFoundAt(tokenI);
            while ((opLvl[tokenI] != -4) || (pl[tokenI] != pLvl))
            {
                if (tokenI == (map.size() - 1))
                {
                    OStringStream description;
                    description << "Can't find closing parenthesis for "
                        << "pow() function.";
                    fatalParseError
                    (
                        index,
                        tl,
                        powFoundAt,
                        tokenI,
                        "equationReader::removePowExponents",
                        description
                    );
                }
                tokenI++;
            }

            const label closeFoundAt(tokenI);
            // Ignore if the exponent is only 1 token
            if ((closeFoundAt - commaFoundAt) > 2)
            {

                // Now create sub-equation
                OStringStream subEqnStream;
                for
                (
                    label subTokenI(commaFoundAt + 1);
                    subTokenI < closeFoundAt;
                    subTokenI++
                )
                {
                    if
                    (
                        tl[subTokenI].isPunctuation()
                     && (tl[subTokenI].pToken() == token::COLON))
                    {
                        subEqnStream << "^";
                    }
                    else
                    {
                        subEqnStream << tl[subTokenI];
                    }
                }
                string subEqnRawText(subEqnStream.str());
                const equation& eqn(operator[](index));
                equation subEqn
                (
                    eqn.name() + "_powExponent_" + name(powFoundAt),
                    subEqnRawText,
                    eqn.overrideDimensions(),
                    eqn.changeDimensions()
                );

                bool eqnCreated(false);
                for (label eqnI(0); eqnI < size(); eqnI++)
                {
                    const equation& eqnTest(operator[](eqnI));
                    if (eqnTest.name() == subEqn.name())
                    {
                        clearEquation(eqnI);
                        eqnTest.setRawText(subEqn.rawText());
                        eqnTest.setOverrideDimensions
                        (
                            subEqn.overrideDimensions()
                        );
                        eqnTest.setChangeDimensions
                        (
                            eqnTest.changeDimensions()
                        );
                        eqnCreated = true;
                    }
                }

                if (!eqnCreated)
                {
                    createEquation(subEqn);
                }

                // Change commaFoundAt + 1 entry to reflect new subEquation
                // reference
                tl[commaFoundAt + 1] = token(subEqn.name());
                map.set
                (
                    commaFoundAt + 1,
                    new equationOperation(findSource(subEqn.name()))
                );
                opLvl[commaFoundAt + 1] = 0;
                pl[commaFoundAt + 1] = pl[commaFoundAt];

                // Remove the subEquation from tl, map, opLvl and pl
                label tokensRemoved(closeFoundAt - (commaFoundAt + 2));
                label newSize(map.size() - tokensRemoved);
                for
                (
                    label subTokenI(commaFoundAt + 2);
                    subTokenI < newSize;
                    subTokenI++
                )
                {
                    tl[subTokenI] = tl[subTokenI + tokensRemoved];
                    map[subTokenI] = map[subTokenI + tokensRemoved];
                    opLvl[subTokenI] = opLvl[subTokenI + tokensRemoved];
                    pl[subTokenI] = pl[subTokenI + tokensRemoved];
                }
                tl.setSize(newSize);
                map.setSize(newSize);
                opLvl.setSize(newSize);
                pl.setSize(newSize);
            }
        }
        tokenI++;
    }
}


Foam::labelList Foam::equationReader::findMaxOperation
(
    const labelList& opLvl,
    const labelList& indices
) const
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
    PtrList<equationOperation>& map,
    labelList& opLvl
) const
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


void Foam::equationReader::findMaxFieldSize(const label equationIndex) const
{
    // *** This function has been made "FULLDEBUG" only since field sizes are
    // *** handled automatically with the evaluateField functions.
    //
    // Look through the operations' sources, find those that are fields or
    // equations, and take their min sizes, by geoIndex.
    // geoIndex is a way to accommodate GeometricFields:
    //      geoIndex = 0: internal field
    //      geoIndex > 0: boundary field
    //          the geoIndex for boundary fields is (hence) 1-indexed.
    //
    // The rules for finding the maximum field size are:
    // 1. sources with a single cell size: scalars, vectors, dictionary, etc.
    //      are assumed to be uniform, and place no limit to the field sizes.
    // 2. The limitations come from field sources and equation sources:
    //      a) The maximum geoIndex is the minimum of these sources.
    //      b) The maximum field size in each geoIndex is the minimum of these
    //          sources.
    // Mismatch rules:
    //  When a mismatch is encountered in geoIndex size or field size, it fails
    //  with a fatalError.  There is one exception that does not fail:
    //      -A regular field and a GeometricField will have mismatching
    //       geoIndices: field = 1, GeometricField > 1.  This is ignored
    //       and we assume the equation is only for the internal field.
#ifdef FULLDEBUG
    labelList maxFieldSizes(0);

    const equation& eqn(operator[](equationIndex));

    forAll(eqn, opindex)
    {
        const equationOperation& op(operator[](equationIndex)[opindex]);
        labelList srcFieldSizes(0);
        switch (op.sourceType())
        {
            case equationOperation::stequation:
            {
                const equation& eqnSrc(operator[](mag(op.sourceIndex() - 1)));
                if (eqnSrc.size() == 0)
                {
                    parse(mag(op.sourceIndex()) - 1);
                }
                srcFieldSizes = eqnSrc.maxFieldSizes();
                break;
            }
            case equationOperation::stscalarFieldSource:
            {
                label sourceIndex(mag(op.sourceIndex() - 1));
                srcFieldSizes.setSize(scalarSources_.geoSize(sourceIndex));
                forAll(srcFieldSizes, geoIndex)
                {
                    srcFieldSizes[geoIndex] = scalarSources_.fieldSize
                    (
                        sourceIndex,
                        geoIndex
                    );
                }
                break;
            }
            case equationOperation::stvectorFieldSource:
            {
                label sourceIndex(mag(op.sourceIndex() - 1));
                srcFieldSizes.setSize(vectorSources_.geoSize(sourceIndex));
                forAll(srcFieldSizes, geoIndex)
                {
                    srcFieldSizes[geoIndex] = vectorSources_.fieldSize
                    (
                        sourceIndex,
                        geoIndex
                    );
                }
                break;
            }
            case equationOperation::sttensorFieldSource:
            {
                label sourceIndex(mag(op.sourceIndex() - 1));
                srcFieldSizes.setSize(tensorSources_.geoSize(sourceIndex));
                forAll(srcFieldSizes, geoIndex)
                {
                    srcFieldSizes[geoIndex] = tensorSources_.fieldSize
                    (
                        sourceIndex,
                        geoIndex
                    );
                }
                break;
            }
            case equationOperation::stdiagTensorFieldSource:
            {
                label sourceIndex(mag(op.sourceIndex() - 1));
                srcFieldSizes.setSize(diagTensorSources_.geoSize(sourceIndex));
                forAll(srcFieldSizes, geoIndex)
                {
                    srcFieldSizes[geoIndex] = diagTensorSources_.fieldSize
                    (
                        sourceIndex,
                        geoIndex
                    );
                }
                break;
            }
            case equationOperation::stsymmTensorFieldSource:
            {
                label sourceIndex(mag(op.sourceIndex() - 1));
                srcFieldSizes.setSize(symmTensorSources_.geoSize(sourceIndex));
                forAll(srcFieldSizes, geoIndex)
                {
                    srcFieldSizes[geoIndex] = symmTensorSources_.fieldSize
                    (
                        sourceIndex,
                        geoIndex
                    );
                }
                break;
            }
            case equationOperation::stsphericalTensorFieldSource:
            {
                label sourceIndex(mag(op.sourceIndex() - 1));
                srcFieldSizes.setSize
                (
                    sphericalTensorSources_.geoSize(sourceIndex)
                );
                forAll(srcFieldSizes, geoIndex)
                {
                    srcFieldSizes[geoIndex] = sphericalTensorSources_.fieldSize
                    (
                        sourceIndex,
                        geoIndex
                    );
                }
                break;
            }
            default:
                // Do nothing
                break;
        }
        if (srcFieldSizes.size() != 0)
        {
            if (maxFieldSizes.size() == 0)
            {
                maxFieldSizes = srcFieldSizes;
                break;
            }
            if (maxFieldSizes.size() == 1)
            {
                if (maxFieldSizes[0] != srcFieldSizes[0])
                {
                    FatalErrorIn("equationReader::findMaxFieldSize")
                        << "Field size mismatch: equation " << eqn.name()
                        << ": fieldSizes [" << maxFieldSizes << "] != "
                        << "fieldSizes [" << srcFieldSizes << "]."
                        << abort(FatalError);
                }
            }
            if
            (
                (maxFieldSizes.size() > 1)
             && (srcFieldSizes.size() == 1)
             && (maxFieldSizes[0] == srcFieldSizes[0])
            )
            {
                // A singleField-only source is encountered.  This equation is
                // now only applicable for the internalField
                maxFieldSizes.setSize(1);
            }
            else if
            (
                (
                    (maxFieldSizes.size() == 1)
                 && (maxFieldSizes[0] != srcFieldSizes[0])
                )
             || (
                    (maxFieldSizes.size() != 1)
                 && (maxFieldSizes != srcFieldSizes)
                )
            )
            {
                FatalErrorIn("equationReader::findMaxFieldSize")
                    << "Field size mismatch: equation " << eqn.name()
                    << ": fieldSizes [" << maxFieldSizes << "] != "
                    << "fieldSizes [" << srcFieldSizes << "]."
                    << abort(FatalError);
            }
        }
    }
    eqn.setMaxFieldSizes(maxFieldSizes);
#endif
}


void Foam::equationReader::dsEqual
(
    dimensionedScalar& dso,
    const dimensionedScalar& dsi
) const
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
)
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
) const
{
    // varName may contain a component: varPart.componentPart
    word varPart(varName);
    word componentPart(word::null);

    // Search for last '.', excluding extreme ends
    for (label i(varName.size() - 2); i > 0; i--)
    {
        if (varName(i, 1) == word("."))
        {
            componentPart = varName(i + 1, varName.size() - (i + 1));
            varPart = varName(0, i);
        }
    }

    // Search order:
    // -other known equations
    // -activeSources;
    // -scalarSources, then scalarFieldSources;
    // -vectorSources, then vectorFieldSources;
    // -tensorSources, then tensorFieldSources;
    // -diagTensorSources, then diagTensorFieldSources;
    // -symmTensorSources, then symmTensorFieldSources;
    // -sphericalTensorSources, then sphericalTensorFieldSources;
    // -dictSources, which may be either:
    //  --a scalar;
    //  --a dimensionedScalar; or
    //  --a yet unknown equation;

    // Searching known equations
    for (label eqs(0); eqs < size(); eqs++)
    {
        if (operator[](eqs).name() == varName)
        {
            return equationOperation
            (
                equationOperation::stequation,
                eqs + 1,
                0,
                0,
                equationOperation::otnone
            );
        }
    }

    forAll(activeSourceNames_, aIndex)
    {
        if (activeSourceNames_[aIndex] == varPart)
        {
            return equationOperation
            (
                equationOperation::stactiveSource,
                aIndex + 1,
                activeSources_[aIndex].lookupComponentIndex(componentPart),
                0,
                equationOperation::otnone
            );
        }
    }

    if (scalarSources_.foundSingle(varPart))
    {
        return equationOperation
        (
            equationOperation::stscalarSource,
            scalarSources_.lookupSingle(varPart) + 1,
            scalarSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (scalarSources_.foundField(varPart))
    {
        return equationOperation
        (
            equationOperation::stscalarFieldSource,
            scalarSources_.lookupField(varPart) + 1,
            scalarSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (vectorSources_.foundSingle(varPart))
    {
        return equationOperation
        (
            equationOperation::stvectorSource,
            vectorSources_.lookupSingle(varPart) + 1,
            vectorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (vectorSources_.foundField(varPart))
    {
        return equationOperation
        (
            equationOperation::stvectorFieldSource,
            vectorSources_.lookupField(varPart) + 1,
            vectorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (tensorSources_.foundSingle(varPart))
    {
        return equationOperation
        (
            equationOperation::sttensorSource,
            tensorSources_.lookupSingle(varPart) + 1,
            tensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (tensorSources_.foundField(varPart))
    {
        return equationOperation
        (
            equationOperation::sttensorFieldSource,
            tensorSources_.lookupField(varPart) + 1,
            tensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (diagTensorSources_.foundSingle(varPart))
    {
        return equationOperation
        (
            equationOperation::stdiagTensorSource,
            diagTensorSources_.lookupSingle(varPart) + 1,
            diagTensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (diagTensorSources_.foundField(varPart))
    {
        return equationOperation
        (
            equationOperation::stdiagTensorFieldSource,
            diagTensorSources_.lookupField(varPart) + 1,
            diagTensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (symmTensorSources_.foundSingle(varPart))
    {
        return equationOperation
        (
            equationOperation::stsymmTensorSource,
            symmTensorSources_.lookupSingle(varPart) + 1,
            symmTensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (symmTensorSources_.foundField(varPart))
    {
        return equationOperation
        (
            equationOperation::stsymmTensorFieldSource,
            symmTensorSources_.lookupField(varPart) + 1,
            symmTensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (sphericalTensorSources_.foundSingle(varPart))
    {
        return equationOperation
        (
            equationOperation::stsphericalTensorSource,
            sphericalTensorSources_.lookupSingle(varPart) + 1,
            sphericalTensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    if (sphericalTensorSources_.foundField(varPart))
    {
        return equationOperation
        (
            equationOperation::stsphericalTensorFieldSource,
            sphericalTensorSources_.lookupField(varPart) + 1,
            sphericalTensorSources_.lookupComponentIndex(componentPart),
            0,
            equationOperation::otnone
        );
    }

    forAll(dictSources_, i)
    {
        if (dictSources_[i].found(varName) && !dictSources_[i].isDict(varName))
        {
            label dictLookupIndex(-1);

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
                    equationOperation::stdictSource,
                    i + 1,
                    0,
                    dictLookupIndex,
                    equationOperation::otnone
                );
            }
            else if (isEquation(srcStrm))
            {
                // Is it a known equation already?
                for (label eqs(0); eqs < size(); eqs++)
                {
                    if (operator[](eqs).name() == varName)
                    {
                        return equationOperation
                        (
                            equationOperation::stequation,
                            eqs + 1,
                            0,
                            0,
                            equationOperation::otnone
                        );
                    }
                }

                // The equation has not been read yet.  Create an unparsed
                // equation.  It will be parsed during evaluate, or if it
                // is later read.
                equation eqn(srcStrm);
                eqn.name() = varName;
                createEquation(eqn);

                return equationOperation
                (
                    equationOperation::stequation,
                    size(),
                    0,
                    0,
                    equationOperation::otnone
                );
            }
        } // end if varName is in dictionary
    } // end dictionary search loop

    // Nothing found, return empty
    return equationOperation
    (
        equationOperation::stnone,
        0,
        0,
        0,
        equationOperation::otnone
    );

}


Foam::label Foam::equationReader::addInternalScalar(const scalar& value) const
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


void Foam::equationReader::setSize(const label newSize) const
{
    eqns_.setSize(newSize);
}

void Foam::equationReader::set(const label index, equation * eqn) const
{
    eqns_.set(index, eqn);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationReader::equationReader(bool showSplash)
:
    scalarSources_
    (
        word("scalar")
    ),
    vectorSources_
    (
        word("vector")
    ),
    tensorSources_
    (
        word("tensor")
    ),
    diagTensorSources_
    (
        word("diagTensor")
    ),
    symmTensorSources_
    (
        word("symmTensor")
    ),
    sphericalTensorSources_
    (
        word("sphericalTensor")
    )
{
    if (showSplash)
    {
        Info
            << "/*                       |---------------------." << token::NL
            << " * This application uses | David L. F. Gaden's |  "
            << "Please cite me if possible" << token::NL
            << " *      .----------------|---------------------'  "
            << "See the wiki for more info" << token::NL
            << " *      | equationReader |  Version:    " << version()
            << token::NL
            << " *      '----------------|  Wiki:       "
            << "github.com/Marupio/equationReader/wiki" << token::NL
            << " */" << endl;
    }
    if (debug)
    {
        reportEmbeddedDispatchFunction_
            = &Foam::equationReader::reportEmbeddedDispatchEnabled;
        reportEmbeddedReturnFunction_
            = &Foam::equationReader::reportEmbeddedReturnEnabled;
    }
    else
    {
        reportEmbeddedDispatchFunction_
            = &Foam::equationReader::reportEmbeddedDispatchDisabled;
        reportEmbeddedReturnFunction_
            = &Foam::equationReader::reportEmbeddedReturnDisabled;
    }
    if ((debug == 1) || (debug == 2) || (debug == 5) || (debug == 6))
    {
        reportScalarEvalStartFunction_
            = &Foam::equationReader::reportScalarEvalStartEnabled;
        reportScalarEvalEndFunction_
            = &Foam::equationReader::reportScalarEvalEndEnabled;
    }
    else
    {
        reportScalarEvalStartFunction_
            = &Foam::equationReader::reportScalarEvalStartDisabled;
        reportScalarEvalEndFunction_
            = &Foam::equationReader::reportScalarEvalEndDisabled;
    }
    if ((debug == 2) || (debug == 6))
    {
        reportScalarOperationFunction_
            = &Foam::equationReader::reportScalarOperationEnabled;
        reportScalarResultFunction_
            = &Foam::equationReader::reportScalarResultEnabled;
    }
    else
    {
        reportScalarOperationFunction_
            = &Foam::equationReader::reportScalarOperationDisabled;
        reportScalarResultFunction_
            = &Foam::equationReader::reportScalarResultDisabled;
    }
    if ((debug == 3) || (debug == 4) || (debug == 5) || (debug == 6))
    {
        reportDimsEvalStartFunction_
            = &Foam::equationReader::reportDimsEvalStartEnabled;
        reportDimsEvalEndFunction_
            = &Foam::equationReader::reportDimsEvalEndEnabled;
    }
    else
    {
        reportDimsEvalStartFunction_
            = &Foam::equationReader::reportDimsEvalStartDisabled;
        reportDimsEvalEndFunction_
            = &Foam::equationReader::reportDimsEvalEndDisabled;
    }
    if ((debug == 4) || (debug == 6))
    {
        reportDimsOperationFunction_
            = &Foam::equationReader::reportDimsOperationEnabled;
        reportDimsResultFunction_
            = &Foam::equationReader::reportDimsResultEnabled;
    }
    else
    {
        reportDimsOperationFunction_
            = &Foam::equationReader::reportDimsOperationDisabled;
        reportDimsResultFunction_
            = &Foam::equationReader::reportDimsResultDisabled;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationReader::~equationReader()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::equationReader::version() const
{
    OStringStream os;
    os << label(equationReaderVersionMajor) << "."
        << label(equationReaderVersionMinor) << "."
        << label(equationReaderVersionBuild);

    return word(os.str());
}



bool Foam::equationReader::found(const word& equationName)
{
    for (label i = 0; i < size(); i++)
    {
        if (operator[](i).name() == equationName)
        {
            return true;
        }
    }
    return false;
}


Foam::label Foam::equationReader::lookup(const word& equationName) const
{
    for (label i = 0; i < size(); i++)
    {
        if (operator[](i).name() == equationName)
        {
            return i;
        }
    }

    return -1;
}


void Foam::equationReader::addSource(const dictionary& dict)
{
    dictSources_.setSize(dictSources_.size() + 1);
    dictSources_.set(dictSources_.size() - 1, &dict);
}


void Foam::equationReader::addSource(const equationVariable& aVar)
{
    label newIndex(activeSources_.size());
    activeSources_.setSize(newIndex + 1);
    activeSources_.set(newIndex, &aVar);
    activeSourceNames_.setSize(newIndex + 1);
    activeSourceNames_[newIndex] = aVar.name();
}


Foam::label Foam::equationReader::readEquation
(
    equation eqn,
    bool okayToReread
)
{
    if (!okayToReread)
    {
        return createEquation(eqn);
    }
    else
    {
        label index(lookup(eqn.name()));
        if (index < 0)
        {
            return createEquation(eqn);
        }
        else if(eqn.rawText() != operator[](index).rawText())
        {
            clearEquation(index);
            return index;
        }
    }
    return -1;
}


Foam::label Foam::equationReader::readEquation
(
    const dictionary& sourceDict,
    const word& eqnName,
    bool okayToReread
)
{
    equation eqn(sourceDict.lookup(eqnName));
    eqn.name() = eqnName;
    return readEquation(eqn);
}


void Foam::equationReader::clearEquation(const label equationIndex) const
{
    operator[](equationIndex).clear();
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
#   ifdef FULLDEBUG
        if ((index < 0) || (index >= size()))
        {
            FatalErrorIn("equationReader::deleteEquation(index)")
                << "Index " << index << " out of bounds (0, " << size() - 1
                << ")"
                << abort(FatalError);
        }
#   endif

    labelList oldToNew(size());
    for (label i(index); i < (size() - 1); i++)
    {
        eqns_[i] = eqns_[i + 1];
    }
    setSize(size() - 1);
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
        is.rewind();
        return true;
    }
    else
    {
        is.rewind();
        return false;
    }
}


bool Foam::equationReader::isWord(ITstream& is)
{
    token firstToken(is);

    if (firstToken.isWord() && is.eof())
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


bool Foam::equationReader::isNamelessDimensionedScalar(ITstream& is)
{
    tokenList tl(10);
    label found(0);
    while (!is.eof())
    {
        tl[found] = token(is);
        found++;
        if (found > 10)
        {
            is.rewind();
            return false;
        }
    }
    if
    (
        (
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
         && tl[9].isNumber()
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
         && tl[7].isNumber()
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


const Foam::equation& Foam::equationReader::operator[]
(
    const label equationIndex
) const
{
    return eqns_[equationIndex];
}


Foam::equation& Foam::equationReader::operator[]
(
    const label equationIndex
)
{
    return eqns_[equationIndex];
}


Foam::label Foam::equationReader::size() const
{
    return eqns_.size();
}

#include "equationReaderAssignFunctionPointers.C"
#include "equationReaderCreateMap.C"
#include "equationReaderEvaluate.C"
#include "equationReaderParse.C"
#include "equationReaderDebugP.C"
#include "equationReaderEvalDimsP.C"
#include "equationReaderEvalScalarP.C"
#include "equationReaderEvalScalarFieldP.C"
#include "equationReaderGetSourceDimsP.C"
#include "equationReaderGetSourceScalarP.C"
#include "equationReaderGetSourceScalarFieldP.C"
// equationReaderIO.C - compile target

// ************************************************************************* //
