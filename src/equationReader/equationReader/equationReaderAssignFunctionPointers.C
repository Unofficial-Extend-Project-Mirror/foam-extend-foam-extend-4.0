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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::equationReader::assignFunctionPointers(const label index) const
{
    forAll(operator[](index), i)
    {
        //Assign getSource functions first
        const equation& eqn(operator[](index));
        equationOperation& eqOp(eqn[i]);
        label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
        switch (operator[](index)[i].sourceType())
        {
            case equationOperation::stnone:
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcNone
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcNone
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcNone
                );
                break;
            case equationOperation::ststorage:
                if (eqOp.operation() == equationOperation::otstore)
                {
                    eqOp.assignSourceDimsFunction
                    (
                        &Foam::equationReader::getDimsSrcNone
                    );
                    eqOp.assignSourceScalarFunction
                    (
                        &Foam::equationReader::getScalarSrcNone
                    );
                    eqOp.assignSourceScalarFieldFunction
                    (
                        &Foam::equationReader::getScalarFieldSrcNone
                    );
                }
                else
                {
                    eqOp.assignSourceDimsFunction
                    (
                        &Foam::equationReader::getDimsSrcStorage
                    );
                    eqOp.assignSourceScalarFunction
                    (
                        &Foam::equationReader::getScalarSrcStorage
                    );
                    eqOp.assignSourceScalarFieldFunction
                    (
                        &Foam::equationReader::getScalarFieldSrcStorage
                    );
                }
                break;
            case equationOperation::stactiveSource:
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcActiveSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcActiveSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcActiveSource
                );
                break;
            case equationOperation::stequation:
                if  (zeroSourceIndex >= size())
                {
                    FatalErrorIn("equationReader::getSouce")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << size() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcEquationCircRefDetect
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcEquationCircRefDetect
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader
                        ::getScalarFieldSrcEquationCircRefDetect
                );
                break;
            case equationOperation::stinternalScalar:
                if (zeroSourceIndex >= internalScalars_.size())
                {
                    FatalErrorIn("equationReader::getSouce")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << internalScalars_.size() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcInternalScalar
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcInternalScalar
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcInternalScalar
                );
                break;
            case equationOperation::stdictSource:
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
                    eqOp.assignSourceDimsFunction
                    (
                        &Foam::equationReader::getDimsSrcDictSourceDScalar
                    );
                    eqOp.assignSourceScalarFunction
                    (
                        &Foam::equationReader::getScalarSrcDictSourceDScalar
                    );
                    eqOp.assignSourceScalarFieldFunction
                    (
                        &Foam::equationReader
                            ::getScalarFieldSrcDictSourceDScalar
                    );
                }
                else if (isScalar(srcStrm))
                {
                    eqOp.assignSourceDimsFunction
                    (
                        &Foam::equationReader::getDimsSrcDictSourceScalar
                    );
                    eqOp.assignSourceScalarFunction
                    (
                        &Foam::equationReader::getScalarSrcDictSourceScalar
                    );
                    eqOp.assignSourceScalarFieldFunction
                    (
                        &Foam::equationReader
                            ::getScalarFieldSrcDictSourceScalar
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
                        << "Expecting a scalar or a dimensionedScalar.  "
                        << "Keyword " << varName << " is referenced by an "
                        << "equation, and therefore can only be one of these "
                        << "two types."
                        << abort(FatalIOError);
                }
                break;
            }
            case equationOperation::stscalarSource:
                if (zeroSourceIndex >= scalarSources_.nSingles())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << scalarSources_.nSingles() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcScalarSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcScalarSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcScalarSource
                );
                break;
            case equationOperation::stscalarFieldSource:
                if (zeroSourceIndex >= scalarSources_.nFields())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << scalarSources_.nFields() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcScalarFieldSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcScalarFieldSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcScalarFieldSource
                );
                break;
            case equationOperation::stvectorSource:
                if (zeroSourceIndex >= vectorSources_.nSingles())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << vectorSources_.nSingles() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcVectorSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcVectorSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcVectorSource
                );
                break;
            case equationOperation::stvectorFieldSource:
                if (zeroSourceIndex >= vectorSources_.nFields())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << scalarSources_.nFields() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcVectorFieldSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcVectorFieldSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcVectorFieldSource
                );
                break;
            case equationOperation::sttensorSource:
                if (zeroSourceIndex >= tensorSources_.nSingles())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << tensorSources_.nSingles() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcTensorSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcTensorSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcTensorSource
                );
                break;
            case equationOperation::sttensorFieldSource:
                if (zeroSourceIndex >= tensorSources_.nFields())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << scalarSources_.nFields() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcTensorFieldSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcTensorFieldSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcTensorFieldSource
                );
                break;
            case equationOperation::stdiagTensorSource:
                if (zeroSourceIndex >= diagTensorSources_.nSingles())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << diagTensorSources_.nSingles() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcDiagTensorSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcDiagTensorSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcDiagTensorSource
                );
                break;
            case equationOperation::stdiagTensorFieldSource:
                if (zeroSourceIndex >= diagTensorSources_.nFields())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << scalarSources_.nFields() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcDiagTensorFieldSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcDiagTensorFieldSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader
                        ::getScalarFieldSrcDiagTensorFieldSource
                );
                break;
            case equationOperation::stsymmTensorSource:
                if (zeroSourceIndex >= symmTensorSources_.nSingles())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << symmTensorSources_.nSingles() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcSymmTensorSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcSymmTensorSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader::getScalarFieldSrcSymmTensorSource
                );
                break;
            case equationOperation::stsymmTensorFieldSource:
                if (zeroSourceIndex >= symmTensorSources_.nFields())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << scalarSources_.nFields() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcSymmTensorFieldSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcSymmTensorFieldSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader
                        ::getScalarFieldSrcSymmTensorFieldSource
                );
                break;
            case equationOperation::stsphericalTensorSource:
                if (zeroSourceIndex >= sphericalTensorSources_.nSingles())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << sphericalTensorSources_.nSingles() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader::getDimsSrcSphericalTensorSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader::getScalarSrcSphericalTensorSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader
                        ::getScalarFieldSrcSphericalTensorSource
                );
                break;
            case equationOperation::stsphericalTensorFieldSource:
                if (zeroSourceIndex >= sphericalTensorSources_.nFields())
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "Index " << zeroSourceIndex << " out of bounds (0, "
                        << scalarSources_.nFields() - 1 << ")"
                        << abort(FatalError);
                }
                eqOp.assignSourceDimsFunction
                (
                    &Foam::equationReader
                        ::getDimsSrcSphericalTensorFieldSource
                );
                eqOp.assignSourceScalarFunction
                (
                    &Foam::equationReader
                        ::getScalarSrcSphericalTensorFieldSource
                );
                eqOp.assignSourceScalarFieldFunction
                (
                    &Foam::equationReader
                        ::getScalarFieldSrcSphericalTensorFieldSource
                );
                break;
        }

        // Assign evaluate functions next
        switch (eqOp.operation())
        {
            case equationOperation::otnone:
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldNone
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarNone
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsNone
                );
                break;
            case equationOperation::otretrieve:
                 eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldRetrieve
                );
                 eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarRetrieve
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsRetrieve
                );
                break;
            case equationOperation::otstore:
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldStore
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarStore
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsStore
                );
                break;
            case equationOperation::otplus:
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldPlus
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarPlus
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsPlusDimCheck
                );
                break;
            case equationOperation::otminus:
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldMinus
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarMinus
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsMinusDimCheck
                );
                break;
            case equationOperation::ottimes:
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldTimes
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarTimes
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsTimes
                );
                break;
            case equationOperation::otdivide:
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldDivide
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarDivide
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsDivide
                );
                break;
            case equationOperation::otpow:
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldPow
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarPow
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsPowDimCheck
                );
                break;
            case equationOperation::otsign:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldSign
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarSign
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsSign
                );
                break;
            case equationOperation::otpos:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldPos
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarPos
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsPos
                );
                break;
            case equationOperation::otneg:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldNeg
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarNeg
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsNeg
                );
                break;
            case equationOperation::otmag:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::assignFunctionPointers")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldMag
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarMag
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsMag
                );
                break;
            case equationOperation::otlimit:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldLimit
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarLimit
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsLimit
                );
                break;
            case equationOperation::otminMod:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldMinMod
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarMinMod
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsMinMod
                );
                break;
            case equationOperation::otsqrtSumSqr:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldSqrtSumSqr
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarSqrtSumSqr
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsSqrtSumSqr
                );
                break;
            case equationOperation::otsqr:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldSqr
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarSqr
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsSqr
                );
                break;
            case equationOperation::otpow3:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldPow3
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarPow3
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsPow3
                );
                break;
            case equationOperation::otpow4:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldPow4
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarPow4
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsPow4
                );
                break;
            case equationOperation::otpow5:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldPow5
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarPow5
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsPow5
                );
                break;
            case equationOperation::otpow6:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldPow6
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarPow6
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsPow6
                );
                break;
            case equationOperation::otinv:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldInv
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarInv
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsInv
                );
                break;
            case equationOperation::otsqrt:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldSqrt
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarSqrt
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsSqrt
                );
                break;
            case equationOperation::otcbrt:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldCbrt
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarCbrt
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsCbrt
                );
                break;
            case equationOperation::othypot:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldHypot
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarHypot
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsHypotDimCheck
                );
                break;
            case equationOperation::otexp:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldExp
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarExp
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsExpDimCheck
                );
                break;
            case equationOperation::otlog:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldLog
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarLog
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsLogDimCheck
                );
                break;
            case equationOperation::otlog10:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldLog10
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarLog10
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsLog10DimCheck
                );
                break;
            case equationOperation::otsin:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldSin
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarSin
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsSinDimCheck
                );
                break;
            case equationOperation::otcos:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldCos
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarCos
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsCosDimCheck
                );
                break;
            case equationOperation::ottan:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldTan
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarTan
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsTanDimCheck
                );
                break;
            case equationOperation::otasin:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldAsin
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarAsin
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsAsinDimCheck
                );
                break;
            case equationOperation::otacos:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldAcos
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarAcos
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsAcosDimCheck
                );
                break;
            case equationOperation::otatan:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldAtan
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarAtan
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsAtanDimCheck
                );
                break;
            case equationOperation::otatan2:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldAtan2
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarAtan2
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsAtan2
                );
                break;
            case equationOperation::otsinh:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldSinh
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarSinh
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsSinhDimCheck
                );
                break;
            case equationOperation::otcosh:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldCosh
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarCosh
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsCoshDimCheck
                );
                break;
            case equationOperation::ottanh:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldTanh
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarTanh
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsTanhDimCheck
                );
                break;
            case equationOperation::otasinh:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldAsinh
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarAsinh
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsAsinhDimCheck
                );
                break;
            case equationOperation::otacosh:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldAcosh
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarAcosh
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsAcoshDimCheck
                );
                break;
            case equationOperation::otatanh:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldAtanh
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarAtanh
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsAtanhDimCheck
                );
                break;
            case equationOperation::oterf:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldErf
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarErf
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsErfDimCheck
                );
                break;
            case equationOperation::oterfc:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldErfc
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarErfc
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsErfcDimCheck
                );
                break;
            case equationOperation::otlgamma:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldLgamma
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarLgamma
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsLgammaDimCheck
                );
                break;
            case equationOperation::otj0:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldJ0
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarJ0
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsJ0DimCheck
                );
                break;
            case equationOperation::otj1:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldJ1
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarJ1
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsJ1DimCheck
                );
                break;
            case equationOperation::otjn:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldJn
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarJn
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsJnDimCheck
                );
                break;
            case equationOperation::oty0:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldY0
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarY0
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsY0DimCheck
                );
                break;
            case equationOperation::oty1:
                if
                (
                    eqOp.sourceType()
                 != equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] takes only one parameter."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldY1
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarY1
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsY1DimCheck
                );
                break;
            case equationOperation::otyn:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldYn
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarYn
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsYnDimCheck
                );
                break;
            case equationOperation::otmax:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldMax
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarMax
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsMaxDimCheck
                );
                break;
            case equationOperation::otmin:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldMin
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarMin
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsMinDimCheck
                );
                break;
            case equationOperation::otstabilise:
                if
                (
                    eqOp.sourceType()
                 == equationOperation::stnone
                )
                {
                    FatalErrorIn("equationReader::evaluate")
                        << "In the equation for "
                        << eqn.name() << ", given by "
                        << token::NL << token::TAB  << eqn.rawText()
                        << token::NL << "Function ["
                        << equationOperation::opName(eqOp.operation())
                        << "] requires two parameters."
                        << abort(FatalError);
                }
                eqOp.assignOpScalarFieldFunction
                (
                    &Foam::equationReader::evalScalarFieldStabilise
                );
                eqOp.assignOpScalarFunction
                (
                    &Foam::equationReader::evalScalarStabilise
                );
                eqOp.assignOpDimsFunction
                (
                    &Foam::equationReader::evalDimsStabilise
                );
                break;
        } // end switch
    } // end forAll equation operations
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
