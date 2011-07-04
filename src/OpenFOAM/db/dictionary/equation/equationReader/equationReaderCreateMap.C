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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::equationReader::createMap
(
    const label index,
    const tokenList& tl,
    equationOperationList& map,
    labelList& opLvl,
    labelList& pl
)
{
//    equation * eqn(&this->operator[](index));
    
    // current parenthesis level - note, a negative parenthesis value indicates
    // that this is the root level of a function, and therefore ',' is allowed
    label p(0);

    forAll(tl, i)
    {
        if (tl[i].isNumber())
        {
            // Internal constant.  Save to internalScalars and record source
            opLvl[i] = 0;
            pl[i] = p;
            map[i] = equationOperation
            (
                equationOperation::slinternalScalar,
                addInternalScalar(tl[i].number()) + 1,
                0,
                equationOperation::otnone
            );
/*
            //Bug fix
            eqn = &this->operator[](index);

            eqn->internalScalars()->setSize
            (
                eqn->internalScalars()->size() + 1
            );
            eqn->internalScalars()->operator[]
            (
                eqn->internalScalars()->size() - 1
            ) = tl[i].number();
            opLvl[i] = 0;
            pl[i] = p;
            map[i] = equationOperation
            (
                equationOperation::slinternalScalar,
                eqn->internalScalars()->size(),
                0,
                equationOperation::otnone
            );
*/
        }
        else if (tl[i].isWord())
        {
            // could be a variable name, function or mathematical constant
            // - check for function first - function is [word][punctuation '(']
            if
            (
                (i < (tl.size() - 1))
             && (tl[i + 1].isPunctuation())
             && (tl[i + 1].pToken() == token::BEGIN_LIST)
            )
            {
                // Function detected; function brackets are negative
                opLvl[i] = 4;
                p = -mag(p) - 1;
                pl[i] = p;
                map[i] = equationOperation
                (
                    equationOperation::slnone,
                    0,
                    0,
                    equationOperation::findOp(tl[i].wordToken())
                );
                    
                if (map[i].operation() == equationOperation::otnone)
                {
                    OStringStream description;
                    description << tl[i].wordToken() << " is not a recognized "
                        << "function.";
                    fatalParseError
                    (
                        index,
                        tl,
                        i,
                        i,
                        "equationReader::parse",
                        description
                    );
                }
                
                // Set next token as well (function opening parenthesis)
                i++;
                opLvl[i] = 4;
                pl[i] = p;
                map[i] = equationOperation
                (
                    equationOperation::slnone,
                    0,
                    0,
                    equationOperation::otnone
                );
            }
            else if
            (
                (tl[i].wordToken() == "e_")
             || (tl[i].wordToken() == "pi_")
             || (tl[i].wordToken() == "twoPi_")
             || (tl[i].wordToken() == "piByTwo_")
             || (tl[i].wordToken() == "GREAT_")
             || (tl[i].wordToken() == "VGREAT_")
             || (tl[i].wordToken() == "ROOTVGREAT_")
             || (tl[i].wordToken() == "SMALL_")
             || (tl[i].wordToken() == "VSMALL_")
             || (tl[i].wordToken() == "ROOTVSMALL_")
            )             
            {
                // Mathematical constant
                
                //Bug fix
//                eqn = &this->operator[](index);

                if
                (
                    findSource(tl[i].wordToken()).sourceList()
                 != equationOperation::slnone
                )
                {
                    // Found a possible conflicting variable name - warn
                    WarningIn("equationReader::createMap")
                        << "Equation for " << eqns_[index].equationName()
                        << ", given by:" << token::NL << token::TAB
                        << eqns_[index].rawText() << token:: NL << "refers "
                        << "to '" << tl[i].wordToken() << "'. Although "
                        << "variable " << tl[i].wordToken() << "was found in "
                        << "the data sources, " << tl[i].wordToken() << " is a"
                        << " mathematical constant. The mathematical constant "
                        << "will be used." << endl;
                }
                
/*                eqn->internalScalars()->setSize
                (
                    eqn->internalScalars()->size() + 1
                );
*/
                opLvl[i] = 0;
                pl[i] = p;
                label internalIndex(0);
                if (tl[i].wordToken() == "e_")
                {
                    internalIndex =
                        addInternalScalar(mathematicalConstant::e) + 1;
                }
                else if (tl[i].wordToken() == "pi_")
                {
                    internalIndex =
                        addInternalScalar(mathematicalConstant::pi) + 1;
                }
                else if (tl[i].wordToken() == "twoPi_")
                {
                    internalIndex =
                        addInternalScalar(mathematicalConstant::twoPi) + 1;
                }
                else if (tl[i].wordToken() == "piByTwo_")
                {
                    internalIndex =
                        addInternalScalar(mathematicalConstant::piByTwo) + 1;
                }
                else if (tl[i].wordToken() == "GREAT_")
                {
                    internalIndex =
                        addInternalScalar(GREAT) + 1;
                }
                else if (tl[i].wordToken() == "VGREAT_")
                {
                    internalIndex =
                        addInternalScalar(VGREAT) + 1;
                }
                else if (tl[i].wordToken() == "ROOTVGREAT_")
                {
                    internalIndex =
                        addInternalScalar(ROOTVGREAT) + 1;
                }
                else if (tl[i].wordToken() == "SMALL_")
                {
                    internalIndex =
                        addInternalScalar(SMALL) + 1;
                }
                else if (tl[i].wordToken() == "VSMALL_")
                {
                    internalIndex =
                        addInternalScalar(VSMALL) + 1;
                }
                else  // tl[i].wordToken() == "ROOTVSMALL_"
                {
                    internalIndex =
                        addInternalScalar(ROOTVSMALL) + 1;
                }
                map[i] = equationOperation
                (
                    equationOperation::slinternalScalar,
                    internalIndex,
                    0,
                    equationOperation::otnone
                );
/*
                if (tl[i].wordToken() == "e_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = mathematicalConstant::e;
                }
                else if (tl[i].wordToken() == "pi_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = mathematicalConstant::pi;
                }
                else if (tl[i].wordToken() == "twoPi_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = mathematicalConstant::twoPi;
                }
                else if (tl[i].wordToken() == "piByTwo_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = mathematicalConstant::piByTwo;
                }
                else if (tl[i].wordToken() == "GREAT_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = GREAT;
                }
                else if (tl[i].wordToken() == "VGREAT_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = VGREAT;
                }
                else if (tl[i].wordToken() == "ROOTVGREAT_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = ROOTVGREAT;
                }
                else if (tl[i].wordToken() == "SMALL_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = SMALL;
                }
                else if (tl[i].wordToken() == "VSMALL_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = VSMALL;
                }
                else if (tl[i].wordToken() == "ROOTVSMALL_")
                {
                    eqn->internalScalars()->operator[]
                    (
                        eqn->internalScalars()->size() - 1
                    ) = ROOTVSMALL;
                }
                opLvl[i] = 0;
                pl[i] = p;
                map[i] = equationOperation
                (
                    equationOperation::slinternalScalar,
                    eqn->internalScalars()->size(),
                    0,
                    equationOperation::otnone
                );
*/
            }
            else
            {
                // Variable name
                opLvl[i] = 0;
                pl[i] = p;
                map[i] = findSource(tl[i].wordToken());
                if (map[i].sourceIndex() == 0)
                {
                    OStringStream description;
                    description << "Variable name " << tl[i].wordToken()
                        << " not found in any available sources.";
                    fatalParseError
                    (
                        index,
                        tl,
                        i,
                        i,
                        "equationReader::parse",
                        description
                    );
                }
            }
        }
        else if (tl[i].isPunctuation())
        {
            switch (tl[i].pToken())
            {
                case token::BEGIN_LIST: // (
                    opLvl[i] = 4;
                    p = mag(p) + 1;
                    pl[i] = p;
                    map[i] = equationOperation
                    (
                        equationOperation::slnone,
                        0,
                        0,
                        equationOperation::otnone
                    );
                    break;
                case token::END_LIST: // )
                {
                    opLvl[i] = -4;
                    pl[i] = p;
                    p = mag(p) - 1;
                    if (p < 0)
                    {
                        OStringStream description;
                        description << "Too many ')'.";
                        fatalParseError
                        (
                            index,
                            tl,
                            i,
                            i,
                            "equationReader::parse",
                            description
                        );

                    }

                    // Look for preceding parenthesis change - was it negative?
                    for (label j(i - 1); j >= 0; j--)
                    {
                        if (mag(pl[j]) == p)
                        {
                            if (pl[j] < 0)
                            {
                                p = -p;
                            }
                            break;
                        }
                    }
                    
                    map[i] = equationOperation
                    (
                        equationOperation::slnone,
                        0,
                        0,
                        equationOperation::otnone
                    );
                    break;
                }
                case token::COMMA: // ,
                    // , is only accepted in a function level parenthesis
                    if (p < 0)
                    {
                        opLvl[i] = 5;
                        pl[i] = p;
                        map[i] = equationOperation
                        (
                            equationOperation::slnone,
                            0,
                            0,
                            equationOperation::otnone
                        );
                    }
                    else
                    {
                        OStringStream description;
                        description << "The comma, ',' does not make sense "
                            << "here.  Only permitted in the root parenthesis "
                            << "level of a function.";
                        fatalParseError
                        (
                            index,
                            tl,
                            i,
                            i,
                            "equationReader::parse",
                            description
                        );
                    }
                    break;
                case token::ADD: // +
                    opLvl[i] = 1;
                    pl[i] = p;
                    map[i] = equationOperation
                    (
                        equationOperation::slnone,
                        0,
                        0,
                        equationOperation::otplus
                    );
                    break;
                case token::SUBTRACT: // -
                    opLvl[i] = 1;
                    pl[i] = p;
                    map[i] = equationOperation
                    (
                        equationOperation::slnone,
                        0,
                        0,
                        equationOperation::otminus
                    );
                    break;
                case token::MULTIPLY: // *
                    opLvl[i] = 2;
                    pl[i] = p;
                    map[i] = equationOperation
                    (
                        equationOperation::slnone,
                        0,
                        0,
                        equationOperation::ottimes
                    );
                    break;
                case token::DIVIDE: // /
                    opLvl[i] = 2;
                    pl[i] = p;
                    map[i] = equationOperation
                    (
                        equationOperation::slnone,
                        0,
                        0,
                        equationOperation::otdivide
                    );
                    break;
                case token::COLON: // :, means ^
                    opLvl[i] = 3;
                    pl[i] = p;
                    map[i] = equationOperation
                    (
                        equationOperation::slnone,
                        0,
                        0,
                        equationOperation::otpow
                    );
                    break;
                default:
                {
                    OStringStream description;
                    description << "Punctuation character '" << tl[i].pToken()
                    << "' is prohibitted.";
                    fatalParseError
                    (
                        index,
                        tl,
                        i,
                        i,
                        "equationReader::parse",
                        description
                    );
                    break;
                }
            } // end punctuation switch
        } // end if punctuation
        else
        {
            OStringStream description;
            description << "Unrecognized token: [" << tl[i] << "].";
            fatalParseError
            (
                index,
                tl,
                i,
                i,
                "equationReader::parse",
                description
            );
        }
    } // mapping loop
    
    if (p)
    {
        OStringStream description;
        description << "Parentheses do not match.  Expecting " << mag(p)
            << " additional ')'s.";
        fatalParseError
        (
            index,
            tl,
            0,
            tl.size() - 1,
            "equationReader::parse",
            description
        );
    }
    
    // Assign negatives (distinguish these from subtraction)
    // The difference is characterized by the preceding character:
    // -preceeded by an operator = negative '+' '-' '*' '/' '^'
    // -preceeded by an open bracket = negative '('
    // -preceeded by a comma = negative ','
    // -preceeded by a variable = subtract 'word' or 'number'
    // -preceeded by a close bracket = subtract ')'
    // Negatives are identified by a negative dictLookupIndex

    if (map[0].operation() == equationOperation::otminus)
    {
        opLvl[0] = 2;
        map[0].dictLookupIndex() = -1;
    }

    for (label i(1); i < map.size(); i++)
    {
        if (map[i].operation() == equationOperation::otminus)
        {
            if (opLvl[i-1] > 0)
            {
               opLvl[i] = 2;
               map[i].dictLookupIndex() = -1;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
