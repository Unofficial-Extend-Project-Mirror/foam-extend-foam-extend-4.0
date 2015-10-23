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

void Foam::equationReader::evalScalarNone
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    FatalErrorIn("equationReader::evalScalarNone")
        << "Empty operation called in equation "
        << operator[](index).name()
        << ", given by:" << token::NL << token::TAB
        << operator[](index).rawText() << token::NL
        << "Empty operations should only exist temporarily during parsing, "
        << "and they should not remain in the operation list at this point.  "
        << "Either you have corrupt data, or this is a bug."
        << abort(FatalError);
}


void Foam::equationReader::evalScalarRetrieve
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = source;
}


void Foam::equationReader::evalScalarStore
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    storeIndex++;
    storageScalars_.setSize(storeIndex + storageOffset + 1);
    storageScalars_[storeIndex + storageOffset] = x;
    x = 0;
}


void Foam::equationReader::evalScalarPlus
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x += source;
}


void Foam::equationReader::evalScalarMinus
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x -= source;
}


void Foam::equationReader::evalScalarTimes
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x *= source;
}


void Foam::equationReader::evalScalarDivide
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x /= source;
}


void Foam::equationReader::evalScalarPow
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = pow(x, source);
}


void Foam::equationReader::evalScalarSign
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = sign(x);
}


void Foam::equationReader::evalScalarPos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = pos(x);
}


void Foam::equationReader::evalScalarNeg
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = neg(x);
}


void Foam::equationReader::evalScalarMag
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = mag(x);
}


void Foam::equationReader::evalScalarLimit
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = limit(x, source);
}


void Foam::equationReader::evalScalarMinMod
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = minMod(x, source);
}


void Foam::equationReader::evalScalarSqrtSumSqr
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = sqrtSumSqr(x, source);
}


void Foam::equationReader::evalScalarSqr
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = sqr(x);
}


void Foam::equationReader::evalScalarPow3
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = pow3(x);
}


void Foam::equationReader::evalScalarPow4
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = pow4(x);
}


void Foam::equationReader::evalScalarPow5
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = pow5(x);
}


void Foam::equationReader::evalScalarPow6
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = pow6(x);
}


void Foam::equationReader::evalScalarInv
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = inv(x);
}


void Foam::equationReader::evalScalarSqrt
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = sqrt(x);
}


void Foam::equationReader::evalScalarCbrt
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = cbrt(x);
}


void Foam::equationReader::evalScalarHypot
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = hypot(x, source);
}


void Foam::equationReader::evalScalarExp
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = exp(x);
}


void Foam::equationReader::evalScalarLog
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = log(x);
}


void Foam::equationReader::evalScalarLog10
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = log10(x);
}


void Foam::equationReader::evalScalarSin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = sin(x);
}

void Foam::equationReader::evalScalarCos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = cos(x);
}


void Foam::equationReader::evalScalarTan
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = tan(x);
}


void Foam::equationReader::evalScalarAsin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = asin(x);
}


void Foam::equationReader::evalScalarAcos
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = acos(x);
}


void Foam::equationReader::evalScalarAtan
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = atan(x);
}


void Foam::equationReader::evalScalarAtan2
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = atan2(x, source);
}


void Foam::equationReader::evalScalarSinh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = sinh(x);
}


void Foam::equationReader::evalScalarCosh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = cosh(x);
}


void Foam::equationReader::evalScalarTanh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = tanh(x);
}


void Foam::equationReader::evalScalarAsinh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = asinh(x);
}


void Foam::equationReader::evalScalarAcosh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = acosh(x);
}


void Foam::equationReader::evalScalarAtanh
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = atanh(x);
}


void Foam::equationReader::evalScalarErf
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = erf(x);
}


void Foam::equationReader::evalScalarErfc
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = erfc(x);
}


void Foam::equationReader::evalScalarLgamma
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = lgamma(x);
}


void Foam::equationReader::evalScalarJ0
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = j0(x);
}


void Foam::equationReader::evalScalarJ1
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = j1(x);
}


void Foam::equationReader::evalScalarJn
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    int xi(x);
    x = jn(xi, source);
}


void Foam::equationReader::evalScalarY0
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = y0(x);
}


void Foam::equationReader::evalScalarY1
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = y1(x);
}


void Foam::equationReader::evalScalarYn
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    int xi(x);
    x = yn(xi, source);
}


void Foam::equationReader::evalScalarMax
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = max(x, source);
}


void Foam::equationReader::evalScalarMin
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = min(x, source);
}


void Foam::equationReader::evalScalarStabilise
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storeIndex,
    scalar& x,
    scalar source
) const
{
    x = stabilise(x, source);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
