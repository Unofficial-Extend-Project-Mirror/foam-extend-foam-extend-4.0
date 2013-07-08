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

void Foam::equationReader::evalNone
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    FatalErrorIn("equationReader::update(index)")
        << "Empty operation called.  Empty operations should only "
        << "exist temporarily during parsing, and they should not "
        << "remain in the operation list at this point.  Either "
        << "you have corrupt data, or this is a bug."
        << abort(FatalError);
}


void Foam::equationReader::evalRetrieve
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, source);
}


void Foam::equationReader::evalRetrieveChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    // Since the last operation of every equation is an otretrieve, this is
    // where we set the dimensions to the overrideDimensions
    ds.name() = source.name();
    ds.value() = source.value();
    ds.dimensions().reset(eqns_[index].overrideDimensions());
}


void Foam::equationReader::evalStore
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    storeIndex++;
    storage_.setSize(storeIndex + storageOffset + 1);
    storage_.set
    (
        storeIndex + storageOffset,
        new dimensionedScalar(ds)
    );
    dsEqual(ds, dimensionedScalar("empty", dimless, 0));
}


void Foam::equationReader::evalPlus
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds = ds + source;
}


void Foam::equationReader::evalPlusDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != source.dimensions())
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalPlusDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    ds = ds + source;
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalPlus
    );
}


void Foam::equationReader::evalPlusChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    source.dimensions().reset(dimless);
    dsEqual(ds, ds + source);
}


void Foam::equationReader::evalMinus
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds = ds - source;
}


void Foam::equationReader::evalMinusDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != source.dimensions())
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalMinusDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    ds = ds - source;
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalMinus
    );
}


void Foam::equationReader::evalMinusChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    source.dimensions().reset(dimless);
    dsEqual(ds, ds - source);
}


void Foam::equationReader::evalTimes
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, ds * source);
}


void Foam::equationReader::evalDivide
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, ds / source);
}


void Foam::equationReader::evalPow
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, pow(ds, source));
}


void Foam::equationReader::evalPowDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalPowDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, pow(ds, source));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalPow
    );
}


void Foam::equationReader::evalPowChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    source.dimensions().reset(dimless);
    dsEqual(ds, pow(ds, source));
}


void Foam::equationReader::evalSign
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds = sign(ds);
}


void Foam::equationReader::evalPos
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds = pos(ds);
}


void Foam::equationReader::evalNeg
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds = neg(ds);
}


void Foam::equationReader::evalMag
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds = mag(ds);
}


void Foam::equationReader::evalLimit
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.value() = limit(ds.value(), source.value());
}


void Foam::equationReader::evalMinMod
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.value() = minMod(ds.value(), source.value());
}


void Foam::equationReader::evalSqrtSumSqr
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.value() = sqrtSumSqr(ds.value(), source.value());
}


void Foam::equationReader::evalSqr
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, sqr(ds));
}


void Foam::equationReader::evalPow3
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, pow3(ds));
}


void Foam::equationReader::evalPow4
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, pow4(ds));
}


void Foam::equationReader::evalPow5
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, pow5(ds));
}


void Foam::equationReader::evalInv
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, pTraits<scalar>::one / ds);
}


void Foam::equationReader::evalSqrt
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, sqrt(ds));
}


void Foam::equationReader::evalCbrt
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, cbrt(ds));
}


void Foam::equationReader::evalHypot
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, hypot(ds, source));
}


void Foam::equationReader::evalHypotDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != source.dimensions())
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalHypotDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, hypot(ds, source));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalHypot
    );
}


void Foam::equationReader::evalHypotChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    source.dimensions().reset(dimless);
    dsEqual(ds, hypot(ds, source));
}


void Foam::equationReader::evalExp
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, exp(ds));
}


void Foam::equationReader::evalExpDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalExpDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, exp(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalExp
    );
}


void Foam::equationReader::evalExpChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, exp(ds));
}


void Foam::equationReader::evalLog
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, log(ds));
}


void Foam::equationReader::evalLogDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalLogDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, log(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalLog
    );
}


void Foam::equationReader::evalLogChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, log(ds));
}


void Foam::equationReader::evalLog10
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, log10(ds));
}


void Foam::equationReader::evalLog10DimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalLog10DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, log10(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalLog10
    );
}


void Foam::equationReader::evalLog10ChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, log10(ds));
}


void Foam::equationReader::evalSin
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, sin(ds));
}


void Foam::equationReader::evalSinDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalSinDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, sin(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalSin
    );
}


void Foam::equationReader::evalSinChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, sin(ds));
}


void Foam::equationReader::evalCos
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, cos(ds));
}


void Foam::equationReader::evalCosDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalCosDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, cos(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalCos
    );
}


void Foam::equationReader::evalCosChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, cos(ds));
}


void Foam::equationReader::evalTan
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, tan(ds));
}


void Foam::equationReader::evalTanDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalTanDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, tan(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalTan
    );
}


void Foam::equationReader::evalTanChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, tan(ds));
}


void Foam::equationReader::evalAsin
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, asin(ds));
}


void Foam::equationReader::evalAsinDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalAsinDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, asin(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalAsin
    );
}


void Foam::equationReader::evalAsinChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, asin(ds));
}


void Foam::equationReader::evalAcos
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, acos(ds));
}


void Foam::equationReader::evalAcosDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalAcosDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, acos(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalAcos
    );
}


void Foam::equationReader::evalAcosChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, acos(ds));
}


void Foam::equationReader::evalAtan
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, atan(ds));
}


void Foam::equationReader::evalAtanDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalAtanDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, atan(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalAtan
    );
}


void Foam::equationReader::evalAtanChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, atan(ds));
}


void Foam::equationReader::evalSinh
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, sinh(ds));
}


void Foam::equationReader::evalSinhDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalSinhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, sinh(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalSinh
    );
}


void Foam::equationReader::evalSinhChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, sinh(ds));
}


void Foam::equationReader::evalCosh
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, cosh(ds));
}


void Foam::equationReader::evalCoshDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalCoshDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, cosh(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalCosh
    );
}


void Foam::equationReader::evalCoshChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, cosh(ds));
}


void Foam::equationReader::evalTanh
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, tanh(ds));
}


void Foam::equationReader::evalTanhDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalTanhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, tanh(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalTanh
    );
}


void Foam::equationReader::evalTanhChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, tanh(ds));
}


void Foam::equationReader::evalAsinh
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, asinh(ds));
}


void Foam::equationReader::evalAsinhDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalAsinhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, asinh(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalAsinh
    );
}


void Foam::equationReader::evalAsinhChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, asinh(ds));
}


void Foam::equationReader::evalAcosh
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, acosh(ds));
}


void Foam::equationReader::evalAcoshDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalAcoshDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, acosh(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalAcosh
    );
}


void Foam::equationReader::evalAcoshChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, acosh(ds));
}


void Foam::equationReader::evalAtanh
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, atanh(ds));
}


void Foam::equationReader::evalAtanhDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalAtanhDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, atanh(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalAtanh
    );
}


void Foam::equationReader::evalAtanhChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, atanh(ds));
}


void Foam::equationReader::evalErf
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, erf(ds));
}


void Foam::equationReader::evalErfDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalErfDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, erf(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalErf
    );
}


void Foam::equationReader::evalErfChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, erf(ds));
}


void Foam::equationReader::evalErfc
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, erfc(ds));
}


void Foam::equationReader::evalErfcDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalErfcDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, erfc(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalErfc
    );
}


void Foam::equationReader::evalErfcChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, erfc(ds));
}


void Foam::equationReader::evalLgamma
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, lgamma(ds));
}


void Foam::equationReader::evalLgammaDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalLgammaDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, lgamma(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalLgamma
    );
}


void Foam::equationReader::evalLgammaChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, lgamma(ds));
}


void Foam::equationReader::evalJ0
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, j0(ds));
}


void Foam::equationReader::evalJ0DimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalJ0DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, j0(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalJ0
    );
}


void Foam::equationReader::evalJ0ChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, j0(ds));
}


void Foam::equationReader::evalJ1
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, j1(ds));
}


void Foam::equationReader::evalJ1DimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalJ1DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, j1(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalJ1
    );
}


void Foam::equationReader::evalJ1ChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, j1(ds));
}


void Foam::equationReader::evalJn
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    int di(ds.value());
    dsEqual(ds, jn(di, source));
}


void Foam::equationReader::evalJnDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (source.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalJnDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    int di(ds.value());
    dsEqual(ds, jn(di, source));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalJn
    );
}


void Foam::equationReader::evalJnChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    int di(ds.value());
    source.dimensions().reset(dimless);
    dsEqual(ds, jn(di, source));
}


void Foam::equationReader::evalY0
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, y0(ds));
}


void Foam::equationReader::evalY0DimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalY0DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, y0(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalY0
    );
}


void Foam::equationReader::evalY0ChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, y0(ds));
}


void Foam::equationReader::evalY1
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, y1(ds));
}


void Foam::equationReader::evalY1DimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (ds.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalY1DimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, y1(ds));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalY1
    );
}


void Foam::equationReader::evalY1ChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    dsEqual(ds, y1(ds));
}


void Foam::equationReader::evalYn
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    int di(ds.value());
    dsEqual(ds, yn(di, source));
}


void Foam::equationReader::evalYnDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        (source.dimensions() != dimless)
     && dimensionSet::debug
    )
    {
        WarningIn("equationReader::evalYnDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    int di(ds.value());
    dsEqual(ds, yn(di, source));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalYn
    );
}


void Foam::equationReader::evalYnChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    source.dimensions().reset(dimless);
    int di(ds.value());
    dsEqual(ds, yn(di, source));
}


void Foam::equationReader::evalMax
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, max(ds, source));
}


void Foam::equationReader::evalMaxDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        dimensionSet::debug
     && (ds.dimensions() != source.dimensions())
    )
    {
        WarningIn("equationReader::evalMaxDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, max(ds, source));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalMax
    );
}


void Foam::equationReader::evalMaxChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.dimensions().reset(dimless);
    source.dimensions().reset(dimless);
    dsEqual(ds, max(ds, source));
}


void Foam::equationReader::evalMin
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    dsEqual(ds, min(ds, source));
}


void Foam::equationReader::evalMinDimCheck
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    if
    (
        dimensionSet::debug
     && (ds.dimensions() != source.dimensions())
    )
    {
        WarningIn("equationReader::evalMinDimCheck")
            << "Dimension error thrown for operation ["
            << equationOperation::opName
            (
                eqns_[index].ops()[i].operation()
            )
            << "] in equation " << eqns_[index].equationName()
            << ", given by:" << token::NL << token::TAB
            << eqns_[index].rawText();
    }
    dsEqual(ds, min(ds, source));
    eqns_[index].ops()[i].assignOpFunction
    (
        &Foam::equationReader::evalMin
    );
}


void Foam::equationReader::evalMinChangeDimensions
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    source.dimensions().reset(dimless);
    ds.dimensions().reset(dimless);
    dsEqual(ds, min(ds, source));
}


void Foam::equationReader::evalStabilise
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    ds.value() = stabilise(ds.value(), source.value());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::equationReader::evaluateField
(
    const label& index,
    scalarList& outputSList,
    dimensionSet& dimensions
)
{
    setListIndex(0);
    dimensionedScalar eval = evaluate(index);
    outputSList[0] = eval.value();
    for(label i(1); i < outputSList.size(); i++)
    {
        setListIndex(i);
        outputSList[i] = evaluate(index).value();
    }
    if
    (
        !eqns_[index].changeDimensions()
     && dimensionSet::debug
     && dimensions != eval.dimensions()
    )
    {
        WarningIn("equationReader::update")
            << "Dimension error thrown for equation "
            << eqns_[index].equationName() << ", given by:"
            << token::NL << token::TAB
            << eqns_[index].rawText();

        dimensions = eval.dimensions();
    }
}


Foam::dimensionedScalar Foam::equationReader::evaluate
(
    const word& equationName
)
{
    label index(lookup(equationName));
    if (index < 0)
    {
        FatalErrorIn("equationReader::evaluate(const word)")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    return evaluate(index);
}


Foam::dimensionedScalar Foam::equationReader::evaluate
(
    const label& index,
    label storageOffset
)
{
    if (debug)
    {
        Info << "Evaluating equation " << eqns_[index].equationName()
            << " at index " << index << ", given by:" << token::NL
            << token::TAB << eqns_[index].rawText() << endl;
    }
#   ifdef FULLDEBUG
    if ((index < 0) || (index >= eqns_.size()))
    {
        FatalErrorIn("equationReader::update(const index)")
            << "Index " << index << " out of bounds (0, " << eqns_.size() - 1
            << ")"
            << abort(FatalError);
    }
#   endif
    
    if (eqns_[index].size() == 0)
    {
        parse(index);
    }

    label storeIndex(-1);
    dimensionedScalar ds("empty", dimless, 0);

    for (label i(0); i < eqns_[index].size(); i++)
    {

#       ifdef FULLDEBUG
        if
        (
            (ds.name() == "empty")
         && (
                eqns_[index].ops()[i].operation()
             != equationOperation::otretrieve
            )
        )
        {
            FatalErrorIn("equationReader::update(index)")
                << "Bad operation list.  Operation at " << i << " either "
                << "follows a 'store', or is the first operation.  Therefore "
                << "it should be retrieve, but it is "
                << eqns_[index].ops()[i].operation() << "."
                << abort(FatalError);
        }
#       endif

        dimensionedScalar source("noSource", dimless, 0);
        
        // Execute getSource function
        dsEqual
        (
            source,
            eqns_[index].ops()[i].getSourceFunction
            (
                this,
                index,
                i,
                storeIndex + storageOffset,
                storageOffset
            )
        );

        // Launch the reportOperationFunction - if debug is greater than one,
        // reportOperationEnabled is called, which posts operation-by-operation
        // information to the console.  Otherwise, reportOperationDiabled is
        // called, which does nothing.
        (*this.*reportOperationFunction_)(index, i, ds);
        
        // Execute the eval function to which this operation points
        eqns_[index].ops()[i].opFunction
        (
            this,
            index,
            i,
            storageOffset,
            storeIndex,
            ds,
            source
        );
        
        // If debug level > 1 this will print the result to the console;
        // otherwise, does nothing.        
        (*this.*reportResultFunction_)(ds);
    }


    ds.name() = eqns_[index].equationName();
    
    //Move one level back up on the dependents_ list
    if (dependents_.size())
    {
        dependents_.setSize(dependents_.size() - 1);
    }    

    storage_.setSize(storageOffset);
    if (debug)
    {
        Info << "Equation evaluated.  Result is: " << ds << endl;
    }
    dsEqual(eqns_[index].lastResult(), ds);
    return ds;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
