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

Foam::dimensionSet Foam::equationReader::evaluateDimsEnabled
(
    const label equationIndex,
    const label maxStoreIndex
) const
{
    return internalEvaluateDimensions(equationIndex, maxStoreIndex);
}


Foam::dimensionSet Foam::equationReader::evaluateDimsDisabled
(
    const label equationIndex,
    const label maxStoreIndex
) const
{
    const equation& eqn(operator[](equationIndex));
    return eqn.overrideDimensions();
}


void Foam::equationReader::evalDimsNone
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    FatalErrorIn("equationReader::evalDimsNone")
        << "Empty operation called in equation "
        << operator[](index).name()
        << ", given by:" << token::NL << token::TAB
        << operator[](index).rawText() << token::NL
        << "Empty operations should only exist temporarily during parsing, "
        << "and they should not remain in the operation list at this point.  "
        << "Either you have corrupt data, or this is a bug."
        << abort(FatalError);
}


void Foam::equationReader::evalDimsRetrieve
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(sourceDims);
}


void Foam::equationReader::evalDimsStore
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    storeIndex++;
    storageDims_.setSize(storeIndex + storageOffset + 1);
    storageDims_.set
    (
        storeIndex + storageOffset,
        new dimensionSet(xDims)
    );
    xDims.reset(dimless);
}


void Foam::equationReader::evalDimsPlus
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    // Do nothing
    //xDims += source;
}


void Foam::equationReader::evalDimsPlusDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        (xDims != sourceDims) && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsPlusDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    xDims += sourceDims;
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsPlus
    );
}


void Foam::equationReader::evalDimsMinus
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims -= sourceDims;
}


void Foam::equationReader::evalDimsMinusDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        (xDims != sourceDims) && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsMinusDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    xDims -= sourceDims;
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsMinus
    );
}


void Foam::equationReader::evalDimsTimes
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims *= sourceDims;
}


void Foam::equationReader::evalDimsDivide
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims /= sourceDims;
}


void Foam::equationReader::evalDimsPow
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    const equation& eqn(operator[](index));
    scalar source
    (
        eqn[i].getSourceScalarFunction
        (
            this,
            index,
            i,
            0,
            0
        )
    );
    xDims = pow(xDims, source);
}


void Foam::equationReader::evalDimsPowDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    const equation& eqn(operator[](index));
    const equationOperation& scalarEqOp(eqn[i]);

    if (scalarEqOp.sourceType() == equationOperation::ststorage)
    {
        FatalErrorIn("equationReader::evalDimsPowCheck")
            << "Bad source for pow() function.  This shouldn't happen and is "
            << "a bug.  Try reducing the exponent part of your pow function "
            << "a single term if possible.  Equation in question is:"
            << operator[](index).name() << ", given by:" << token::NL
            << token::TAB << operator[](index).rawText()
            << abort(FatalError);
    }

    if
    (
        !sourceDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsPowDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    scalar source
    (
        eqn[i].getSourceScalarFunction
        (
            this,
            index,
            i,
            0,
            0
        )
    );
    xDims = pow(xDims, source);
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsPow
    );
}


void Foam::equationReader::evalDimsSign
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(sign(xDims));
}


void Foam::equationReader::evalDimsPos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(pos(xDims));
}


void Foam::equationReader::evalDimsNeg
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(neg(xDims));
}


void Foam::equationReader::evalDimsMag
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(mag(xDims));
}


void Foam::equationReader::evalDimsLimit
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    // Do nothing
}


void Foam::equationReader::evalDimsMinMod
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    // Do nothing
}


void Foam::equationReader::evalDimsSqrtSumSqr
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    // Do nothing
}


void Foam::equationReader::evalDimsSqr
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(sqr(xDims));
}


void Foam::equationReader::evalDimsPow3
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(pow3(xDims));
}


void Foam::equationReader::evalDimsPow4
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(pow4(xDims));
}


void Foam::equationReader::evalDimsPow5
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(pow5(xDims));
}


void Foam::equationReader::evalDimsPow6
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(pow6(xDims));
}


void Foam::equationReader::evalDimsInv
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(inv(xDims));
}


void Foam::equationReader::evalDimsSqrt
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(sqrt(xDims));
}


void Foam::equationReader::evalDimsCbrt
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(cbrt(ds).dimensions());
}


void Foam::equationReader::evalDimsHypot
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(hypot(xDims, sourceDims));
}


void Foam::equationReader::evalDimsHypotDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        (xDims != sourceDims) && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsHypotDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 1.0);
    dimensionedScalar dsSource("tempSource", sourceDims, 1.0);
    xDims.reset(hypot(ds, dsSource).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsHypot
    );
}


void Foam::equationReader::evalDimsExp
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(exp(xDims));
}


void Foam::equationReader::evalDimsExpDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsExpDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(exp(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsExp
    );
}


void Foam::equationReader::evalDimsLog
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(log(xDims));
}


void Foam::equationReader::evalDimsLogDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsLogDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(log(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsLog
    );
}


void Foam::equationReader::evalDimsLog10
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(log10(xDims));
}


void Foam::equationReader::evalDimsLog10DimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsLog10DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(log10(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsLog10
    );
}


void Foam::equationReader::evalDimsSin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsSinDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsSinDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(sin(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsSin
    );
}


void Foam::equationReader::evalDimsCos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsCosDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsCosDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(cos(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsCos
    );
}


void Foam::equationReader::evalDimsTan
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsTanDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsTanDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(tan(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsTan
    );
}


void Foam::equationReader::evalDimsAsin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsAsinDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsAsinDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(asin(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsAsin
    );
}


void Foam::equationReader::evalDimsAcos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsAcosDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsAcosDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(acos(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsAcos
    );
}


void Foam::equationReader::evalDimsAtan
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsAtanDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsAtanDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(atan(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsAtan
    );
}


void Foam::equationReader::evalDimsAtan2
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    xDims.reset(dimless);
}


void Foam::equationReader::evalDimsSinh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsSinhDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsSinhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(sinh(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsSinh
    );
}


void Foam::equationReader::evalDimsCosh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsCoshDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsCoshDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(cosh(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsCosh
    );
}


void Foam::equationReader::evalDimsTanh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsTanhDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsTanhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(tanh(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsTanh
    );
}


void Foam::equationReader::evalDimsAsinh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsAsinhDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsAsinhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(asinh(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsAsinh
    );
}


void Foam::equationReader::evalDimsAcosh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsAcoshDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsAcoshDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 1.0);
    xDims.reset(acosh(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsAcosh
    );
}


void Foam::equationReader::evalDimsAtanh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsAtanhDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsAtanhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(atanh(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsAtanh
    );
}


void Foam::equationReader::evalDimsErf
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsErfDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsErfDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(erf(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsErf
    );
}


void Foam::equationReader::evalDimsErfc
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsErfcDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsErfcDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(erfc(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsErfc
    );
}


void Foam::equationReader::evalDimsLgamma
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsLgammaDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsLgammaDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 1.0);
    xDims.reset(lgamma(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsLgamma
    );
}


void Foam::equationReader::evalDimsJ0
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsJ0DimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsJ0DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(j0(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsJ0
    );
}


void Foam::equationReader::evalDimsJ1
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsJ1DimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsJ1DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 0.0);
    xDims.reset(j1(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsJ1
    );
}


void Foam::equationReader::evalDimsJn
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsJnDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        (!sourceDims.dimensionless() || !xDims.dimensionless())
     && dimensionSet::debug
    )
    {
        FatalErrorIn("equationReader::evalDimsJnDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText() << token::NL
            << "jn(a, b) - a and b must be dimensionless."
            << abort(FatalError);
    }
    xDims.reset(dimless);
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsJn
    );
}


void Foam::equationReader::evalDimsY0
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsY0DimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsY0DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 1.0);
    xDims.reset(y0(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsY0
    );
}


void Foam::equationReader::evalDimsY1
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsY1DimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        !xDims.dimensionless() && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalDimsY1DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    dimensionedScalar ds("temp", xDims, 1.0);
    xDims.reset(y1(ds).dimensions());
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsY1
    );
}


void Foam::equationReader::evalDimsYn
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(dimless);
}


void Foam::equationReader::evalDimsYnDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        (!sourceDims.dimensionless() || !xDims.dimensionless())
     && dimensionSet::debug
    )
    {
        FatalErrorIn("equationReader::evalDimsYnDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText() << token::NL
            << "yn(a, b) - a and b must be dimensionless."
            << abort(FatalError);
    }
    xDims.reset(dimless);
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsYn
    );
}


void Foam::equationReader::evalDimsMax
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(max(xDims, sourceDims));
}


void Foam::equationReader::evalDimsMaxDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        dimensionSet::debug
     && (xDims != sourceDims)
    )
    {
        WarningIn("equationReader::evalMaxDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    xDims.reset(max(xDims, sourceDims));
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsMax
    );
}


void Foam::equationReader::evalDimsMin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    //Do nothing
    //xDims.reset(min(xDims, sourceDims));
}


void Foam::equationReader::evalDimsMinDimCheck
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    if
    (
        dimensionSet::debug
     && (xDims != sourceDims)
    )
    {
        WarningIn("equationReader::evalMinDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                operator[](index)[i].operation()
            )
            << "] in equation " << operator[](index).name()
            << ", given by:" << token::NL << token::TAB
            << operator[](index).rawText();
    }
    xDims.reset(min(xDims, sourceDims));
    operator[](index)[i].assignOpDimsFunction
    (
        &Foam::equationReader::evalDimsMin
    );
}


void Foam::equationReader::evalDimsStabilise
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    // Do nothing
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
