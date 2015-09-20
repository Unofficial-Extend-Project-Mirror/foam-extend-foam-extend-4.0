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

void Foam::equationReader::evalScalarFieldNone
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    FatalErrorIn("equationReader::evalScalarFieldNone")
        << "Empty operation called in equation "
        << operator[](index).name()
        << ", given by:" << token::NL << token::TAB
        << operator[](index).rawText() << token::NL
        << "Empty operations should only exist temporarily during parsing, "
        << "and they should not remain in the operation list at this point.  "
        << "Either you have corrupt data, or this is a bug."
        << abort(FatalError);
}


void Foam::equationReader::evalScalarFieldRetrieve
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    x = source;
}


void Foam::equationReader::evalScalarFieldStore
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    storeIndex++;
    storageScalarFields_.setSize(storeIndex + storageOffset + 1);
    storageScalarFields_.set
    (
        storeIndex + storageOffset,
        new scalarField(x)
    );
    x = 0.0;
}


void Foam::equationReader::evalScalarFieldPlus
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    x += source;
}


void Foam::equationReader::evalScalarFieldMinus
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    x -= source;
}


void Foam::equationReader::evalScalarFieldTimes
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    x *= source;
}


void Foam::equationReader::evalScalarFieldDivide
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    x /= source;
}


void Foam::equationReader::evalScalarFieldPow
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    pow(x, x, source);
}


void Foam::equationReader::evalScalarFieldSign
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    sign(x, x);
}


void Foam::equationReader::evalScalarFieldPos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    pos(x, x);
}


void Foam::equationReader::evalScalarFieldNeg
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    neg(x, x);
}


void Foam::equationReader::evalScalarFieldMag
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    x = mag(x);
}


void Foam::equationReader::evalScalarFieldLimit
(
    const equationReader * eqnReader,
    const label index,
    const label iIn,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    TFOR_ALL_F_OP_FUNC_F_F
    (
        scalar, x, =, ::Foam::limit, scalar, x, scalar, source
    )
}


void Foam::equationReader::evalScalarFieldMinMod
(
    const equationReader * eqnReader,
    const label index,
    const label iIn,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    TFOR_ALL_F_OP_FUNC_F_F
    (
        scalar, x, =, ::Foam::minMod, scalar, x, scalar, source
    )
}


void Foam::equationReader::evalScalarFieldSqrtSumSqr
(
    const equationReader * eqnReader,
    const label index,
    const label iIn,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    TFOR_ALL_F_OP_FUNC_F_F
    (
        scalar, x, =, ::Foam::sqrtSumSqr, scalar, x, scalar, source
    )
}


void Foam::equationReader::evalScalarFieldSqr
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    sqr(x, x);
}


void Foam::equationReader::evalScalarFieldPow3
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    pow3(x, x);
}


void Foam::equationReader::evalScalarFieldPow4
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    pow4(x, x);
}


void Foam::equationReader::evalScalarFieldPow5
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    pow5(x, x);
}


void Foam::equationReader::evalScalarFieldPow6
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    pow6(x, x);
}


void Foam::equationReader::evalScalarFieldInv
(
    const equationReader * eqnReader,
    const label index,
    const label iIn,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    TFOR_ALL_F_OP_FUNC_F
    (
        scalar, x, =, ::Foam::inv, scalar, x
    )
}


void Foam::equationReader::evalScalarFieldSqrt
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    sqrt(x, x);
}


void Foam::equationReader::evalScalarFieldCbrt
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    pow(x, x, 1.0/3.0);
}


void Foam::equationReader::evalScalarFieldHypot
(
    const equationReader * eqnReader,
    const label index,
    const label iIn,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    TFOR_ALL_F_OP_FUNC_F_F
    (
        scalar, x, =, ::Foam::hypot, scalar, x, scalar, source
    )
}


void Foam::equationReader::evalScalarFieldExp
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    exp(x, x);
}


void Foam::equationReader::evalScalarFieldLog
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    log(x, x);
}


void Foam::equationReader::evalScalarFieldLog10
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    log10(x, x);
}


void Foam::equationReader::evalScalarFieldSin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    sin(x, x);
}

void Foam::equationReader::evalScalarFieldCos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    cos(x, x);
}


void Foam::equationReader::evalScalarFieldTan
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    tan(x, x);
}


void Foam::equationReader::evalScalarFieldAsin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    asin(x, x);
}


void Foam::equationReader::evalScalarFieldAcos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    acos(x, x);
}


void Foam::equationReader::evalScalarFieldAtan
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    atan(x, x);
}


void Foam::equationReader::evalScalarFieldAtan2
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    atan2(x, x, source);
}


void Foam::equationReader::evalScalarFieldSinh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    sinh(x, x);
}


void Foam::equationReader::evalScalarFieldCosh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    cosh(x, x);
}


void Foam::equationReader::evalScalarFieldTanh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    tanh(x, x);
}


void Foam::equationReader::evalScalarFieldAsinh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    asinh(x, x);
}


void Foam::equationReader::evalScalarFieldAcosh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    acosh(x, x);
}


void Foam::equationReader::evalScalarFieldAtanh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    atanh(x, x);
}


void Foam::equationReader::evalScalarFieldErf
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    erf(x, x);
}


void Foam::equationReader::evalScalarFieldErfc
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    erfc(x, x);
}


void Foam::equationReader::evalScalarFieldLgamma
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    lgamma(x, x);
}


void Foam::equationReader::evalScalarFieldJ0
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    j0(x, x);
}


void Foam::equationReader::evalScalarFieldJ1
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    j1(x, x);
}


void Foam::equationReader::evalScalarFieldJn
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    int xi(x[0]);
    jn(x, xi, source);
}


void Foam::equationReader::evalScalarFieldY0
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    y0(x, x);
}


void Foam::equationReader::evalScalarFieldY1
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    y1(x, x);
}


void Foam::equationReader::evalScalarFieldYn
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    int xi(x[0]);
    yn(x, xi, source);
}


void Foam::equationReader::evalScalarFieldMax
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    max(x, x, source);
}


void Foam::equationReader::evalScalarFieldMin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    min(x, x, source);
}


void Foam::equationReader::evalScalarFieldStabilise
(
    const equationReader * eqnReader,
    const label index,
    const label iIn,
    const label storageOffset,
    label& storeIndex,
    scalarField& x,
    const scalarField& source
) const
{
    TFOR_ALL_F_OP_FUNC_F_F
    (
        scalar, x, =, ::Foam::stabilise, scalar, x, scalar, source
    )
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
