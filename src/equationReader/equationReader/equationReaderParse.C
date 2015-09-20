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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::equationReader::parse(label index) const
{
    if (debug)
    {
        Info << "Parsing equation " << operator[](index).name()
            << " at index " << index << "." << endl;
    }
    if ((index > size()) || (index < 0))
    {
        FatalErrorIn("equationReader::parse(index)")
            << "Index " << index << " out of bounds (0, "
            << size() - 1 << ")"
            << abort(FatalError);
    }

    clearEquation(index);

    // First, ensure there are no ':' or '&' characters.  This is to
    // accommodate the stringPreconditioner, which uses these special
    // characters as work-arounds to limitations imposed by the token class.
    if
    (
        operator[](index).rawText().string::find(":") != string::npos
     && operator[](index).rawText().string::find("&") != string::npos
    )
    {
        FatalErrorIn("equationReader::parse")
            << "Parsing error in the equation for "
            << operator[](index).name() << ", given by:" << token::NL
            << token::NL << token::TAB << operator[](index).rawText()
            << token::NL << token::NL << "Colons, ':', and ampersands '&' are "
            << "prohibitted."
            << abort(FatalError);
    }

    // Precondition the string, and load it into a stream
    IStringStream rawStream(stringPreconditioner(operator[](index).rawText()));
    tokenList tl;

    bool forceEnd(false);

    // Read tokens from raw equation string stream
    while (!rawStream.eof() && !forceEnd)
    {
        tl.setSize(tl.size() + 1);
        tl[tl.size() - 1] = token(rawStream);

        // Bug fix - equations ending in brackets read an extra token of type
        // FATALERROR at the end, caused by string replace ')' with ' ) ' above
        if (tl[tl.size() - 1].type() == token::FATALERROR)
        {
            tl.setSize(tl.size() - 1);
            forceEnd = true;
        }
    }

    // map:
    // - variable / constant: conatins source data only (first three fields)
    // - operation: conatains operation number only (last field)
    // - brackets, comma: all are zero
    PtrList<equationOperation> map(tl.size());

    // opLvl: level of operation precedence
    // 0. variable
    // 1. + -
    // 2. * / and negatives
    // 3. ^
    // 4. ( ) and functions; special case -4 indicates close bracket ')'
    // 5. , used only in functions
    labelList opLvl(tl.size());

    // parenthesis level, negative means function root
    labelList pl(tl.size());

    createMap(index, tl, map, opLvl, pl);

    // As a bug fix, we find any instance of pow(a,b) and remove the b part to
    // a seperate equation.  Necessary because we perform scalar and dimension
    // evaluations seperately. pow is the only function where we can't figure
    // out the dimensions without the value of the exponent.  The only way we
    // are guaranteed to know the scalar value is if it is a stand-alone
    // equation.
    removePowExponents(index, tl, map, opLvl, pl);

    /* This is useful for debugging, so I left it in
    if (debug)
    {
        Info << "tokenList: " << endl;
        forAll(tl, i)
        {
            Info << tl[i];
            if (tl[i].isNumber())
            {
                Info << " isNumber ";
            }
            if (tl[i].isPunctuation())
            {
                Info << " isPunctuation ";
            }
            if (tl[i].isWord())
            {
                Info << " isWord ";
            }
            Info << endl;
        }
        Info << "opLvl is: " << opLvl << endl;
        Info << "pl is: " << pl << endl;
        Info << "Map is: " << map << endl;
    }*/

    // In the main parsing loop, we create labelLists of indices that specify
    // what part of the equation we are working with.  As the equation is
    // parsed, the indices lists will shrink, but tl, pl, opLvl, and map will
    // not.  Indices lists:
    //
    // - eqnIndices - the full working list
    // - subEqnIndices - the current group with the highest parenthesis level
    // - subEqnIndices2 - used when evaluating multiparameter functions
    //
    // Anytime we are 'trimming', we are removing elements from these lists.
    //
    // The main parsing loop:
    // - find the max parenthesis level magnitude
    // - if pl < 0, it is a function, look for a comma
    // -- no comma:
    //    send the expression to parseExpression
    //    if it is a function, evaluate the function
    //    store the result, trim the indices
    // -- comma
    //    same as above, except parseExpression both sides of the comma
    // - once the eqnIndices are down to a single size, parsing is done

    label storeIndex(-1);

    // Create an index list of all the tokens we're working with - initially
    // it is all of them
    labelList eqnIndices(tl.size());
    forAll(eqnIndices, i)
    {
        eqnIndices[i] = i;
    }

    // Main parsing loop
    while (eqnIndices.size() > 1)
    {
        labelList subEqnIndices(findMaxParenthesis(pl, eqnIndices));
        if (opLvl[subEqnIndices[0]] == 4)
        {
            // Expression is enclosed in brackets - trim them
            if (pl[subEqnIndices[0]] < 0)
            {
                // This is a function:
                // - trim function name, but leave it in parent indexList
                // - first bracket is second index
                pl[subEqnIndices[0]] = 0;
                trimList(subEqnIndices, 0, 0);
                trimListWithParent(eqnIndices, subEqnIndices, 0, 0);
            }
            else
            {
                // Not a function, first bracket is first index
                trimListWithParent(eqnIndices, subEqnIndices, 0, 0);
            }

            // Trimming trailing bracket
            trimListWithParent
            (
                eqnIndices,
                subEqnIndices,
                subEqnIndices.size() - 1,
                subEqnIndices.size() - 1
            );
        }

        // Move negatives into the source index
        absorbNegatives(index, tl, eqnIndices, subEqnIndices, map, opLvl);

        label commaIndex(-1);
        label commaPos(-1);

        // Look for a comma
        forAll(subEqnIndices, i)
        {
            if (opLvl[subEqnIndices[i]] == 5)
            {
                commaIndex = i;
                commaPos = subEqnIndices[i];
                break;
            }
        }
        if (subEqnIndices.size() == 2)
        {
            OStringStream description;
            description << "Empty expression '()' found.";
            fatalParseError
            (
                index,
                tl,
                subEqnIndices[0],
                subEqnIndices[subEqnIndices.size() - 1],
                "equationReader::parse",
                description
            );
        }

        if (commaIndex == -1)
        {
            // standard parenthesis or single parameter function
            label resultIndex
            (
                parseExpression
                (
                    index,
                    tl,
                    opLvl,
                    map,
                    subEqnIndices,
                    storeIndex
                )
            );

            trimListWithParent
            (
                eqnIndices,
                subEqnIndices,
                0,
                subEqnIndices.size() - 1,
                findIndex(resultIndex, subEqnIndices)
            );

            label currentIndex(-1);

            if (pl[resultIndex] < 0)
            {
                // This is a single parameter function call - evaluate it
                operator[](index).setSize(operator[](index).size() + 3);

                // retrieve parameter value
                operator[](index).set
                (
                    operator[](index).size() - 3,
                    new equationOperation
                    (
                        map[resultIndex].sourceType(),
                        map[resultIndex].sourceIndex(),
                        map[resultIndex].componentIndex(),
                        map[resultIndex].dictLookupIndex(),
                        equationOperation::otretrieve
                    )
                );

                // perform function operation
                currentIndex = findIndex(resultIndex, eqnIndices);
                operator[](index).set
                (
                    operator[](index).size() - 2,
                    new equationOperation
                    (
                        equationOperation::stnone,
                        0,
                        0,
                        0,
                        map[eqnIndices[currentIndex - 1]].operation()
                    )
                );

                // store result
                storeIndex++;
                operator[](index).set
                (
                    operator[](index).size() - 1,
                    new equationOperation
                    (
                        equationOperation::ststorage,
                        storeIndex + 1,
                        0,
                        0,
                        equationOperation::otstore
                    )
                );

                // Set term in map to store location
                map[resultIndex] =
                    operator[](index)
                    [
                        operator[](index).size() - 1
                    ];
                map[resultIndex].operation() = equationOperation::otnone;

                // Trim function call from indices list
                currentIndex = findIndex(resultIndex, eqnIndices);
                trimList(eqnIndices, currentIndex - 1, currentIndex - 1);
            }
            // reduce the parenthesis level of result
            pl[resultIndex] = mag(pl[resultIndex]) - 1;

            if (pl[resultIndex] > 0)
            {
                // Look for preceding parenthesis change - was it negative?
                currentIndex = findIndex(resultIndex, eqnIndices);
                for (label i(currentIndex - 1); i >= 0; i--)
                {
                    if (mag(pl[eqnIndices[i]]) == pl[eqnIndices[currentIndex]])
                    {
                        if (pl[eqnIndices[i]] < 0)
                        {
                            pl[eqnIndices[currentIndex]] =
                                -pl[eqnIndices[currentIndex]];
                        }
                        break;
                    }
                }
            }
        } // end standard parenthesis / single parameter function
        else if
        (
            (commaIndex < 1) || (commaIndex >= (subEqnIndices.size() - 1))
        )
        {
            OStringStream description;
            description << "Misplaced comma.  '(,[expression)' or "
                << "'([expression],)' found.";
            fatalParseError
            (
                index,
                tl,
                commaPos,
                commaPos,
                "equationReader::parse",
                description
            );
        }
        else
        {
            // multi-parameter function
            // Split the expression into two - before & after the comma
            labelList subEqnIndices2
            (
                subEqnIndices.size() - commaIndex - 1
            );
            forAll(subEqnIndices2, i)
            {
                subEqnIndices2[i] = subEqnIndices[i + commaIndex + 1];
            }
            subEqnIndices.setSize(commaIndex + 1);
            trimListWithParent
            (
                eqnIndices,
                subEqnIndices,
                commaIndex,
                commaIndex
            );

            // Parse the first parameter
            label resultIndex
            (
                parseExpression
                (
                    index,
                    tl,
                    opLvl,
                    map,
                    subEqnIndices,
                    storeIndex
                )
            );

            trimListWithParent
            (
                eqnIndices,
                subEqnIndices,
                0,
                subEqnIndices.size() - 1,
                findIndex(resultIndex, subEqnIndices)
            );

            // Parse the second parameter
            label resultIndex2
            (
                parseExpression
                (
                    index,
                    tl,
                    opLvl,
                    map,
                    subEqnIndices2,
                    storeIndex
                )
            );

            trimListWithParent
            (
                eqnIndices,
                subEqnIndices2,
                0,
                subEqnIndices2.size() - 1,
                findIndex(resultIndex2, subEqnIndices2)
            );

            // Perform multiparameter function operations
            // first retrieve the first parameter
            operator[](index).setSize(operator[](index).size() + 3);
            operator[](index).set
            (
                operator[](index).size() - 3,
                new equationOperation
                (
                        map[resultIndex].sourceType(),
                        map[resultIndex].sourceIndex(),
                        map[resultIndex].componentIndex(),
                        map[resultIndex].dictLookupIndex(),
                        equationOperation::otretrieve
                )
            );

            // perform the function operation (2nd parameter is source)
            label currentIndex(findIndex(resultIndex, eqnIndices));
            operator[](index).set
            (
                operator[](index).size() - 2,
                new equationOperation
                (
                    map[resultIndex2].sourceType(),
                    map[resultIndex2].sourceIndex(),
                    map[resultIndex2].componentIndex(),
                    map[resultIndex2].dictLookupIndex(),
                    map[eqnIndices[currentIndex - 1]].operation()
                )
            );

            // store result
            storeIndex++;
            operator[](index).set
            (
                operator[](index).size() - 1,
                new equationOperation
                (
                    equationOperation::ststorage,
                    storeIndex + 1,
                    0,
                    0,
                    equationOperation::otstore
                )
            );

            // Set term in map to store location
            map.set
            (
                resultIndex,
                new equationOperation
                (
                    operator[](index)[operator[](index).size() - 1]
                )
            );
            map[resultIndex].operation() = equationOperation::otnone;

            // trim function call from indices list
            trimList(eqnIndices, currentIndex - 1, currentIndex - 1);

            // trim second parameter from indices list
            label currentIndex2(findIndex(resultIndex2, eqnIndices));
            trimList(eqnIndices, currentIndex2, currentIndex2);

            // reduce the parenthesis level of result
            pl[resultIndex] = mag(pl[resultIndex]) - 1;

            if (pl[resultIndex] > 0)
            {
                currentIndex = findIndex(resultIndex, eqnIndices);
                // Look for preceding parenthesis change - was it negative?
                for (label i(currentIndex - 1); i >= 0; i--)
                {
                    if (mag(pl[eqnIndices[i]]) == pl[eqnIndices[currentIndex]])
                    {
                        if (pl[eqnIndices[i]] < 0)
                        {
                            pl[eqnIndices[currentIndex]] =
                                -pl[eqnIndices[currentIndex]];
                        }
                        break;
                    }
                }
            }
            // break;
        } // end default case
    } // end main parsing loop

    // Special case - equations with only one effective term:
    // e.g. "2", "-2", "-(2)", "-(((((2)))))", etc..
    // will complete their parse with an empty operation list
    if (operator[](index).size() == 0)
    {

        labelList subEqnIndices(findMaxParenthesis(pl, eqnIndices));
        absorbNegatives(index, tl, eqnIndices, subEqnIndices, map, opLvl);

        if (opLvl[subEqnIndices[0]] != 0)
        {
            OStringStream description;
            description << "Expecting a variable or literal constant.";
            fatalParseError
            (
                index,
                tl,
                0,
                0,
                "equationReader::parse",
                description
            );
        }

        // Add two operations - the first retrieves the variable, the second
        // is a dummy because the last operation is trimmed before exitting
        // equationReader::parse
        operator[](index).setSize(operator[](index).size() + 2);

        // retrieve parameter value
        operator[](index).set
        (
            operator[](index).size() - 2,
            new equationOperation
            (
                map[subEqnIndices[0]].sourceType(),
                map[subEqnIndices[0]].sourceIndex(),
                map[subEqnIndices[0]].componentIndex(),
                map[subEqnIndices[0]].dictLookupIndex(),
                equationOperation::otretrieve
            )
        );

        // Store this result
        operator[](index).set
        (
            operator[](index).size() - 1,
            new equationOperation
            (
                equationOperation::ststorage,
                storeIndex + 1,
                0,
                0,
                equationOperation::otstore
            )
        );
    }

    // The last operation is an otstore.  Add an otretrieve to finalize.
    // We could eliminate the last otstore, but this will miss the final
    // absorbNegatives, if one existed.
    operator[](index).setSize(operator[](index).size() + 1);
    operator[](index).set
    (
        operator[](index).size() - 1,
        new equationOperation
        (
            map[eqnIndices[0]].sourceType(),
            map[eqnIndices[0]].sourceIndex(),
            map[eqnIndices[0]].componentIndex(),
            map[eqnIndices[0]].dictLookupIndex(),
            equationOperation::otretrieve
        )
    );

    // Link the eval function pointers (for efficiency)
    assignFunctionPointers(index);

    // Determine the maximum field sizes available for this equation
    findMaxFieldSize(index);

    if (debug)
    {
        Info << "Parsing complete for equation " << operator[](index).name()
            << " at index " << index << "." << endl;
    }
}


Foam::label Foam::equationReader::parseExpression
(
    label index,
    const tokenList& tl,
    const labelList& opLvl,
    PtrList<equationOperation>& map,
    const labelList& indices,
    label& storeIndex
) const
{
//    equation * eqn(&this->operator[](index));
    label returnMe(-1);
    labelList subIndices(indices);
    labelList opIndices;

    bool done(false);
    while (!done)
    {
        opIndices = findMaxOperation(opLvl, subIndices);
        if (!(opIndices.size() % 2))
        {
            OStringStream description;
            description << "Expected pattern: [number] [operation] [number] "
                << "[operation] ... violated.";
            fatalParseError
            (
                index,
                tl,
                opIndices[0],
                opIndices[opIndices.size() - 1],
                "equationReader::parseExpression",
                description
            );
        }
        if (!opIndices.size())
        {
            OStringStream description;
            description << "Empty expression found, e.g. '()'.";
            fatalParseError
            (
                index,
                tl,
                subIndices[0],
                subIndices[subIndices.size() - 1],
                "equationReader::parseExpression",
                description
            );
        }
        if (opIndices.size() == 1)
        {
            // This means only one term existed between the brackets, nothing
            // needs to be done.
            if (opLvl[opIndices[0]] != 0)
            {
                OStringStream description;
                description << "Detected an isolated operator, e.g. (+).";
                fatalParseError
                (
                    index,
                    tl,
                    opIndices[0],
                    opIndices[0],
                    "equationReader::parse",
                    description
                );
            }
            done = true;
            returnMe = opIndices[0];
        }
        else if (opIndices.size() > 1)
        {
            // More than one term.  Do a retrieve, then enter the operations
            // loop.
            if (opLvl[opIndices[0]] != 0)
            {
                OStringStream description;
                description << "Expected pattern: [number] [operation] "
                    << "[number] [operation] ... violated.  First token is "
                    << "not a [number].";
                fatalParseError
                (
                    index,
                    tl,
                    opIndices[0],
                    opIndices[opIndices.size() - 1],
                    "equationReader::parseExpression",
                    description
                );
            }
            operator[](index).setSize(operator[](index).size() + 1);

            operator[](index).set
            (
                operator[](index).size() - 1,
                new equationOperation
                (
                    map[opIndices[0]].sourceType(),
                    map[opIndices[0]].sourceIndex(),
                    map[opIndices[0]].componentIndex(),
                    map[opIndices[0]].dictLookupIndex(),
                    equationOperation::otretrieve
                )
            );

            trimListWithParent(subIndices, opIndices, 0, 0);

            // Begin operations loop
            while (opIndices.size() > 1)
            {
                operator[](index).setSize(operator[](index).size() + 1);
                operator[](index).set
                (
                    operator[](index).size() - 1,
                    new equationOperation
                    (
                        map[opIndices[1]].sourceType(),
                        map[opIndices[1]].sourceIndex(),
                        map[opIndices[1]].componentIndex(),
                        map[opIndices[1]].dictLookupIndex(),
                        map[opIndices[0]].operation()
                    )
                );

                if (opIndices.size() > 2)
                {
                    // Remove the two operation indices from the working list
                    trimListWithParent(subIndices, opIndices, 0, 1);
                }
                else
                {
                    // Last operation, perform a store
                    operator[](index).setSize(operator[](index).size() + 1);
                    storeIndex++;
                    operator[](index).set
                    (
                        operator[](index).size() - 1,
                        new equationOperation
                        (
                            equationOperation::ststorage,
                            storeIndex + 1,
                            0,
                            0,
                            equationOperation::otstore
                        )
                    );

                    returnMe = opIndices[1];
                    map.set
                    (
                        opIndices[1],
                        new equationOperation
                        (
                            operator[](index)[operator[](index).size() - 1]
                        )
                    );
                    map[opIndices[1]].operation() = equationOperation::otnone;
                    trimListWithParent(subIndices, opIndices, 0, 0);
                }
            } // end operations loop
        } // end if (opIndices.size() > 1)
    } // main parsing loop
    return returnMe;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
