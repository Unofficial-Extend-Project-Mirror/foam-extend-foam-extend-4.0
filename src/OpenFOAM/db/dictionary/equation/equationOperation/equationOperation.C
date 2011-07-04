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

#include "dimensionedScalar.H"
#include "equationReader.H"
#include "equationOperation.H"
//#include "equationOperationList.H"

class dimensionedScalar;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::equationOperation::typeName = "equationOperation";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationOperation::equationOperation()
{}


Foam::equationOperation::equationOperation
(
    sourceListType sourceList,
    label sourceIndex,
    label dictLookupIndex,
    operationType operation,
    dimensionedScalar (Foam::equationReader::*getSourceFunction)
    (
        equationReader *,
        const label,
        const label,
        const label,
        const label
    ),
    void (Foam::equationReader::*opFunction)
    (
        equationReader *,
        const label,
        const label,
        const label,
        label&,
        dimensionedScalar&,
        dimensionedScalar
    )
)
:
    sourceList_(sourceList),
    sourceIndex_(sourceIndex),
    dictLookupIndex_(dictLookupIndex),
    operation_(operation),
    getSourceFunction_(getSourceFunction),
    opFunction_(opFunction)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationOperation::~equationOperation()
{}

// * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * //
Foam::equationOperation::operationType
    Foam::equationOperation::findOp(const Foam::word& opName)
{
    if (opName == "retrieve")
    {
        return otretrieve;
    }
    else if (opName == "store")
    {
        return otstore;
    }
    else if (opName == "plus")
    {
        return otplus;
    }
    else if (opName == "minus")
    {
        return otminus;
    }
    else if (opName == "times")
    {
        return ottimes;
    }
    else if (opName == "divide")
    {
        return otdivide;
    }
    else if (opName == "pow")
    {
        return otpow;
    }
    else if (opName == "sign")
    {
        return otsign;
    }
    else if (opName == "pos")
    {
        return otpos;
    }
    else if (opName == "neg")
    {
        return otneg;
    }
    else if (opName == "mag")
    {
        return otmag;
    }
    else if (opName == "limit")
    {
        return otlimit;
    }
    else if (opName == "minMod")
    {
        return otminMod;
    }
    else if (opName == "sqrtSumSqr")
    {
        return otsqrtSumSqr;
    }
    else if (opName == "sqr")
    {
        return otsqr;
    }
    else if (opName == "pow3")
    {
        return otpow3;
    }
    else if (opName == "pow4")
    {
        return otpow4;
    }
    else if (opName == "pow5")
    {
        return otpow5;
    }
    else if (opName == "inv")
    {
        return otinv;
    }
    else if (opName == "sqrt")
    {
        return otsqrt;
    }
    else if (opName == "cbrt")
    {
        return otcbrt;
    }
    else if (opName == "hypot")
    {
        return othypot;
    }
    else if (opName == "exp")
    {
        return otexp;
    }
    else if (opName == "log")
    {
        return otlog;
    }
    else if (opName == "log10")
    {
        return otlog10;
    }
    else if (opName == "sin")
    {
        return otsin;
    }
    else if (opName == "cos")
    {
        return otcos;
    }
    else if (opName == "tan")
    {
        return ottan;
    }
    else if (opName == "asin")
    {
        return otasin;
    }
    else if (opName == "acos")
    {
        return otacos;
    }
    else if (opName == "atan")
    {
        return otatan;
    }
    else if (opName == "sinh")
    {
        return otsinh;
    }
    else if (opName == "cosh")
    {
        return otcosh;
    }
    else if (opName == "tanh")
    {
        return ottanh;
    }
    else if (opName == "asinh")
    {
        return otasinh;
    }
    else if (opName == "acosh")
    {
        return otacosh;
    }
    else if (opName == "atanh")
    {
        return otatanh;
    }
    else if (opName == "erf")
    {
        return oterf;
    }
    else if (opName == "erfc")
    {
        return oterfc;
    }
    else if (opName == "lgamma")
    {
        return otlgamma;
    }
    else if (opName == "j0")
    {
        return otj0;
    }
    else if (opName == "j1")
    {
        return otj1;
    }
    else if (opName == "jn")
    {
        return otjn;
    }
    else if (opName == "y0")
    {
        return oty0;
    }
    else if (opName == "y1")
    {
        return oty1;
    }
    else if (opName == "yn")
    {
        return otyn;
    }
    else if (opName == "max")
    {
        return otmax;
    }
    else if (opName == "min")
    {
        return otmin;
    }
    else if (opName == "stabilise")
    {
        return otstabilise;
    }
    else
    {
        return otnone;
    }
}


Foam::word Foam::equationOperation::opName
(
    const Foam::equationOperation::operationType& op
)
{
    switch (op)
    {
        case otnone:
            return "none";
        case otretrieve:
            return "retrieve";
        case otstore:
            return "store";
        case otplus:
            return "plus";
        case otminus:
            return "minus";
        case ottimes:
            return "times";
        case otdivide:
            return "divide";
        case otpow:
            return "pow";
        case otsign:
            return "sign";
        case otpos:
            return "pos";
        case otneg:
            return "neg";
        case otmag:
            return "mag";
        case otlimit:
            return "limit";
        case otminMod:
            return "minMod";
        case otsqrtSumSqr:
            return "sqrtSumSqr";
        case otsqr:
            return "sqr";
        case otpow3:
            return "pow3";
        case otpow4:
            return "pow4";
        case otpow5:
            return "pow5";
        case otinv:
            return "inv";
        case otsqrt:
            return "sqrt";
        case otcbrt:
            return "cbrt";
        case othypot:
            return "hypot";
        case otexp:
            return "exp";
        case otlog:
            return "log";
        case otlog10:
            return "log10";
        case otsin:
            return "sin";
        case otcos:
            return "cos";
        case ottan:
            return "tan";
        case otasin:
            return "asin";
        case otacos:
            return "acos";
        case otatan:
            return "atan";
        case otsinh:
            return "sinh";
        case otcosh:
            return "cosh";
        case ottanh:
            return "tanh";
        case otasinh:
            return "asinh";
        case otacosh:
            return "acosh";
        case otatanh:
            return "atanh";
        case oterf:
            return "erf";
        case oterfc:
            return "erfc";
        case otlgamma:
            return "lgamma";
        case otj0:
            return "j0";
        case otj1:
            return "j1";
        case otjn:
            return "jn";
        case oty0:
            return "y0";
        case oty1:
            return "y1";
        case otyn:
            return "yn";
        case otmax:
            return "max";
        case otmin:
            return "min";
        case otstabilise:
            return "stabilise";
    }
}


Foam::word Foam::equationOperation::sourceName
(
    const Foam::equationOperation::sourceListType& sl
)
{
    switch (sl)
    {
        case slnone:
            return "none";
        case sldictSource:
            return "dictionary";
        case slexternalDScalar:
            return "dimensionedScalar";
        case slexternalScalar:
            return "scalar";
        case slexternalScalarList:
            return "scalarList";
        case slinternalScalar:
            return "constant";
        case slequation:
            return "equation";
        case slstorage:
            return "memory";
    }
}


void Foam::equationOperation::assignSourceFunction
(
    dimensionedScalar (Foam::equationReader::*getSourceFunction)
    (
        equationReader *,
        const label,
        const label,
        const label,
        const label
    )
)
{
    getSourceFunction_ = getSourceFunction;
}


void Foam::equationOperation::assignOpFunction
(
    void (Foam::equationReader::*opFunction)
    (
        equationReader *,
        const label,
        const label,
        const label,
        label&,
        dimensionedScalar&,
        dimensionedScalar
    )
)
{
    opFunction_ = opFunction;
}


Foam::dimensionedScalar Foam::equationOperation::getSourceFunction
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    return (eqnReader->*getSourceFunction_)
    (
        eqnReader,
        equationIndex,
        equationOperationIndex,
        maxStoreIndex,
        storageOffset
    );
}


void Foam::equationOperation::opFunction
(
    equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storageIndex,
    dimensionedScalar& ds,
    dimensionedScalar source
)
{
    (eqnReader->*opFunction_)
    (
        eqnReader,
        index,
        i,
        storageOffset,
        storageIndex,
        ds,
        source
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

int Foam::operator==(const equationOperation& I1, const equationOperation& I2)
{
    return
    (
        I1.sourceList() == I2.sourceList()
     && I1.sourceIndex() == I2.sourceIndex()
     && I1.dictLookupIndex() == I2.dictLookupIndex()
     && I1.operation() == I2.operation()
    );
}


int Foam::operator!=(const equationOperation& I1, const equationOperation& I2)
{
    // Invert the '==' operator ('0'='false')
    return I1 == I2 ? 0 : 1;
}


// * * * * * * * * * * * * * Friend IOstream Operators * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, equationOperation& I)
{
    label sl(I.sourceList_);
    label op(I.operation_);
    
    is >> sl >> I.sourceIndex_ >> I.dictLookupIndex_ >> op;
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const equationOperation& I)
{
    label sl(I.sourceList_);
    label op(I.operation_);

    return os   << nl << "/* sourceList: */\t" << sl << nl
                << "/* sourceIndex:*/\t" << I.sourceIndex_ << nl
                << "/* dictIndex   */\t" << I. dictLookupIndex_ << nl
                << "/* operation:  */\t" << op << nl;
}


// ************************************************************************* //
