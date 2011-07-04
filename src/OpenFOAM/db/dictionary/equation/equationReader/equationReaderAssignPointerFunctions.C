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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::equationReader::assignFunctionPointers(const label index)
{
    forAll(eqns_[index], i)
    {
        //Assign getSource functions first
        equation& eqn(eqns_[index]);
        equationOperation& eqOp(eqn[i]);
        label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
        switch (eqns_[index].ops()[i].sourceList())
        {
            case equationOperation::slnone:
                eqOp.assignSourceFunction
                (
                    &Foam::equationReader::getSourceNone
                );
                break;
            case equationOperation::sldictSource:
            {
                if (zeroSourceIndex >= dictSources_.size())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << dictSources_.size() - 1 << ")"
                        << abort(FatalError);
                }
                
                word varName(dictLookups_[eqOp.dictLookupIndex()]);
                
                ITstream srcStrm
                (
                    dictSources_[zeroSourceIndex].lookup(varName)
                );
                if (isDimensionedScalar(srcStrm))
                {
                    eqOp.assignSourceFunction
                    (
                        &Foam::equationReader::getSourceDictSourceDScalar
                    );
                }
                else if (isScalar(srcStrm))
                {
                    eqOp.assignSourceFunction
                    (
                        &Foam::equationReader::getSourceDictSourceScalar
                    );
                }
                else
                {
                    // Neither scalar nor dimensionedScalar
                    FatalIOErrorIn
                    (
                        "equationReader::assignFunctionPointers",
                        dictSources_[zeroSourceIndex]
                    )
                        << "Expecting a scalar or a dimensionedScalar.  Keyword "
                        << varName << " is referenced by an equation, and therfore"
                        << " can only be one of these two types."
                        << exit(FatalIOError);
                }
                break;
            }
            case equationOperation::slexternalDScalar:
                if (zeroSourceIndex >= externalDScalars_.size())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << externalDScalars_.size() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceFunction
                (
                    &Foam::equationReader::getSourceExternalDScalar
                );
                break;
            case equationOperation::slexternalScalar:
                if (zeroSourceIndex >= externalScalars_.size())
                {
                    FatalErrorIn("equationReader::assignPointerFunctions")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << externalScalars_.size() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceFunction
                (
                    &Foam::equationReader::getSourceExternalScalar
                );
                break;
            case equationOperation::slexternalScalarList:
                if (zeroSourceIndex >= externalScalarLists_.size())
                {
                    FatalErrorIn("equationReader::assignPointerFunctions")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << externalScalarLists_.size() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceFunction
                (
                    &Foam::equationReader::getSourceExternalScalarList
                );
                break;
            case equationOperation::slinternalScalar:
                if (zeroSourceIndex >= internalScalars_.size())
                {
                    FatalErrorIn("equationReader::getSouce")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << internalScalars_.size() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceFunction
                (
                    &Foam::equationReader::getSourceInternalScalar
                );
                break;
            case equationOperation::slequation:
                if  (zeroSourceIndex >= eqns_.size())
                {
                    FatalErrorIn("equationReader::getSouce")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << eqns_.size() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceFunction
                (
                    &Foam::equationReader::getSourceEquationCircRefDetect
                );
                break;
            case equationOperation::slstorage:
                if (eqOp.operation() == equationOperation::otstore)
                {
                    eqOp.assignSourceFunction
                    (
                        &Foam::equationReader::getSourceNone
                    );
                }
                else
                {
                    eqOp.assignSourceFunction
                    (
                        &Foam::equationReader::getSourceStorage
                    );
                }
                break;
        }

        // Assign evaluate functions next
        switch (eqOp.operation())
        {
            case equationOperation::otnone:
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalNone
                );
                break;
            case equationOperation::otretrieve:
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalRetrieveChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalRetrieve
                    );
                }
                break;
            case equationOperation::otstore:
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalStore
                );
                break;
            case equationOperation::otplus:
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalPlusChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalPlusDimCheck
                    );
                }
                break;
            case equationOperation::otminus:
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalMinusChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalMinusDimCheck
                    );
                }
                break;
            case equationOperation::ottimes:
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalTimes
                );
                break;
            case equationOperation::otdivide:
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalDivide
                );
                break;
            case equationOperation::otpow:
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalPowChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalPowDimCheck
                    );
                }
                break;
            case equationOperation::otsign:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'sign' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalSign
                );
                break;
            case equationOperation::otpos:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'pos' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalPos
                );
                break;
            case equationOperation::otneg:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'neg' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalNeg
                );
                break;
            case equationOperation::otmag:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'mag' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalMag
                );
                break;
            case equationOperation::otlimit:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'limit' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalLimit
                );
                break;
            case equationOperation::otminMod:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'minMod' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalMinMod
                );
                break;
            case equationOperation::otsqrtSumSqr:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'sqrtSumSqr' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalSqrtSumSqr
                );
                break;
            case equationOperation::otsqr:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'sqr' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalSqr
                );
                break;
            case equationOperation::otpow3:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'pow3' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalPow3
                );
                break;
            case equationOperation::otpow4:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'pow4' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalPow4
                );
                break;
            case equationOperation::otpow5:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'pow5' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalPow5
                );
                break;
            case equationOperation::otinv:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'inv' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalInv
                );
                break;
            case equationOperation::otsqrt:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'sqrt' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalSqrt
                );
                break;
            case equationOperation::otcbrt:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'cbrt' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalCbrt
                );
                break;
            case equationOperation::othypot:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'hypot' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalHypotChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalHypotDimCheck
                    );
                }
                break;
            case equationOperation::otexp:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'exp' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalExpChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalExpDimCheck
                    );
                }
                break;
            case equationOperation::otlog:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'log' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalLogChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalLogDimCheck
                    );
                }
                break;
            case equationOperation::otlog10:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'log10' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalLog10ChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalLog10DimCheck
                    );
                }
                break;
            case equationOperation::otsin:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'sin' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalSinChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalSinDimCheck
                    );
                }
                break;
            case equationOperation::otcos:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'cos' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalCosChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalCosDimCheck
                    );
                }
                break;
            case equationOperation::ottan:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'tan' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalTanChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalTanDimCheck
                    );
                }
                break;
            case equationOperation::otasin:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'asin' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAsinChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAsinDimCheck
                    );
                }
                break;
            case equationOperation::otacos:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'acos' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAcosChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAcosDimCheck
                    );
                }
                break;
            case equationOperation::otatan:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'atan' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAtanChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAtanDimCheck
                    );
                }
                break;
            case equationOperation::otsinh:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'sinh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalSinhChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalSinhDimCheck
                    );
                }
                break;
            case equationOperation::otcosh:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'cosh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalCoshChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalCoshDimCheck
                    );
                }
                break;
            case equationOperation::ottanh:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'tanh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalTanhChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalTanhDimCheck
                    );
                }
                break;
            case equationOperation::otasinh:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'sinh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAsinhChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAsinhDimCheck
                    );
                }
                break;
            case equationOperation::otacosh:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'acosh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAcoshChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAcoshDimCheck
                    );
                }
                break;
            case equationOperation::otatanh:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'atanh' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAtanhChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalAtanhDimCheck
                    );
                }
                break;
            case equationOperation::oterf:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'erf' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalErfChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalErfDimCheck
                    );
                }
                break;
            case equationOperation::oterfc:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'erfc' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalErfcChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalErfcDimCheck
                    );
                }
                break;
            case equationOperation::otlgamma:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'lgamma' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalLgammaChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalLgammaDimCheck
                    );
                }
                break;
            case equationOperation::otj0:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'j0' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalJ0ChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalJ0DimCheck
                    );
                }
                break;
            case equationOperation::otj1:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'j1' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalJ1ChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalJ1DimCheck
                    );
                }
                break;
            case equationOperation::otjn:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'cbrt' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalJnChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalJnDimCheck
                    );
                }
                break;
            case equationOperation::oty0:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'y0' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalY0ChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalY0DimCheck
                    );
                }
                break;
            case equationOperation::oty1:
                if
                (
                    eqOp.sourceList()
                 != equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'y1' takes only one "
                        << "parameter."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalY1ChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalY1DimCheck
                    );
                }
                break;
            case equationOperation::otyn:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'yn' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalY1ChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalY1DimCheck
                    );
                }
                break;
            case equationOperation::otmax:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'max' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalMaxChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalMaxDimCheck
                    );
                }
                break;
            case equationOperation::otmin:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'min' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                if (eqn.changeDimensions())
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalMinChangeDimensions
                    );
                }
                else
                {
                    eqOp.assignOpFunction
                    (
                        &Foam::equationReader::evalMinDimCheck
                    );
                }
                break;
            case equationOperation::otstabilise:
                if
                (
                    eqOp.sourceList()
                 == equationOperation::slnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.equationName() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function 'stabilise' requires two "
                        << "parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpFunction
                (
                    &Foam::equationReader::evalStabilise
                );
                break;
        } // end switch
    } // end forAll equation operations
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
