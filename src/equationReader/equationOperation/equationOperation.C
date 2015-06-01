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


Foam::equationOperation::equationOperation(const equationOperation& eqop)
:
    source_(eqop.source_),
    sourceIndex_(eqop.sourceIndex_),
    componentIndex_(eqop.componentIndex_),
    dictLookupIndex_(eqop.dictLookupIndex_),
    operation_(eqop.operation_),
    getSourceScalarFieldFunction_(eqop.getSourceScalarFieldFunction_),
    opScalarFieldFunction_(eqop.opScalarFieldFunction_),
    getSourceScalarFunction_(eqop.getSourceScalarFunction_),
    opScalarFunction_(eqop.opScalarFunction_),
    getSourceDimsFunction_(eqop.getSourceDimsFunction_),
    opDimsFunction_(eqop.opDimsFunction_)
{}


Foam::equationOperation::equationOperation
(
    sourceTypeEnum source,
    label sourceIndex,
    label componentIndex,
    label dictLookupIndex,
    operationType operation,
    const scalarField& (Foam::equationReader::*getSourceScalarFieldFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        const label
    ) const,
    void (Foam::equationReader::*opScalarFieldFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        label&,
        scalarField&,
        const scalarField&
    ) const,
    scalar (Foam::equationReader::*getSourceScalarFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        const label
    ) const,
    void (Foam::equationReader::*opScalarFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        label&,
        scalar&,
        scalar
    ) const,
    dimensionSet (Foam::equationReader::*getSourceDimsFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        const label
    ) const,
    void (Foam::equationReader::*opDimsFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        label&,
        dimensionSet&,
        dimensionSet
    ) const
)
:
    source_(source),
    sourceIndex_(sourceIndex),
    componentIndex_(componentIndex),
    dictLookupIndex_(dictLookupIndex),
    operation_(operation),
    getSourceScalarFieldFunction_(getSourceScalarFieldFunction),
    opScalarFieldFunction_(opScalarFieldFunction),
    getSourceScalarFunction_(getSourceScalarFunction),
    opScalarFunction_(opScalarFunction),
    getSourceDimsFunction_(getSourceDimsFunction),
    opDimsFunction_(opDimsFunction)
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
    else if (opName == "pow6")
    {
        return otpow6;
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
    else if (opName == "atan2")
    {
        return otatan2;
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
        case otpow6:
            return "pow6";
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
        case otatan2:
            return "atan2";
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
        default:
            return "unlisted";
    }
}


Foam::word Foam::equationOperation::sourceName
(
    const Foam::equationOperation::sourceTypeEnum& st
)
{
    switch (st)
    {
        case stnone:
            return "none";
        case ststorage:
            return "memory";
        case stactiveSource:
            return "activeEquationVariable";
        case stequation:
            return "equation";
        case stinternalScalar:
            return "constant";
        case stdictSource:
            return "dictionary";
        case stscalarSource:
            return "scalar";
        case stscalarFieldSource:
            return "scalarField";
        case stvectorSource:
            return "vector";
        case stvectorFieldSource:
            return "vectorField";
        case sttensorSource:
            return "tensor";
        case sttensorFieldSource:
            return "tensorField";
        case stdiagTensorSource:
            return "diagTensor";
        case stdiagTensorFieldSource:
            return "diagTensorField";
        case stsymmTensorSource:
            return "symmTensor";
        case stsymmTensorFieldSource:
            return "symmTensorField";
        case stsphericalTensorSource:
            return "sphericalTensor";
        case stsphericalTensorFieldSource:
            return "sphericalTensorField";
        default:
            return "unlisted";
    }
}


void Foam::equationOperation::assignSourceScalarFieldFunction
(
    const scalarField& (Foam::equationReader::*getSourceScalarFieldFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        const label
    ) const
) const
{
    getSourceScalarFieldFunction_ = getSourceScalarFieldFunction;
}


void Foam::equationOperation::assignOpScalarFieldFunction
(
    void (Foam::equationReader::*opScalarFieldFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        label&,
        scalarField&,
        const scalarField&
    ) const
) const
{
    opScalarFieldFunction_ = opScalarFieldFunction;
}


void Foam::equationOperation::assignSourceScalarFunction
(
    scalar (Foam::equationReader::*getSourceScalarFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        const label
    ) const
) const
{
    getSourceScalarFunction_ = getSourceScalarFunction;
}


void Foam::equationOperation::assignOpScalarFunction
(
    void (Foam::equationReader::*opScalarFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        label&,
        scalar&,
        scalar
    ) const
) const
{
    opScalarFunction_ = opScalarFunction;
}


void Foam::equationOperation::assignSourceDimsFunction
(
    dimensionSet (Foam::equationReader::*getSourceDimsFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        const label
    ) const
) const
{
    getSourceDimsFunction_ = getSourceDimsFunction;
}


void Foam::equationOperation::assignOpDimsFunction
(
    void (Foam::equationReader::*opDimsFunction)
    (
        const equationReader *,
        const label,
        const label,
        const label,
        label&,
        dimensionSet&,
        dimensionSet
    ) const
) const
{
    opDimsFunction_ = opDimsFunction;
}


const Foam::scalarField& Foam::equationOperation::getSourceScalarFieldFunction
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    return (eqnReader->*getSourceScalarFieldFunction_)
    (
        eqnReader,
        equationIndex,
        equationOperationIndex,
        maxStoreIndex,
        storageOffset
    );
}


void Foam::equationOperation::opScalarFieldFunction
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storageIndex,
    scalarField& x,
    const scalarField& source
) const
{
    (eqnReader->*opScalarFieldFunction_)
    (
        eqnReader,
        index,
        i,
        storageOffset,
        storageIndex,
        x,
        source
    );
}


Foam::scalar Foam::equationOperation::getSourceScalarFunction
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    return (eqnReader->*getSourceScalarFunction_)
    (
        eqnReader,
        equationIndex,
        equationOperationIndex,
        maxStoreIndex,
        storageOffset
    );
}


void Foam::equationOperation::opScalarFunction
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storageIndex,
    scalar& x,
    scalar source
) const
{
    (eqnReader->*opScalarFunction_)
    (
        eqnReader,
        index,
        i,
        storageOffset,
        storageIndex,
        x,
        source
    );
}


Foam::dimensionSet Foam::equationOperation::getSourceDimsFunction
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    return (eqnReader->*getSourceDimsFunction_)
    (
        eqnReader,
        equationIndex,
        equationOperationIndex,
        maxStoreIndex,
        storageOffset
    );
}


void Foam::equationOperation::opDimsFunction
(
    const equationReader * eqnReader,
    const label index,
    const label i,
    const label storageOffset,
    label& storageIndex,
    dimensionSet& xDims,
    dimensionSet sourceDims
) const
{
    (eqnReader->*opDimsFunction_)
    (
        eqnReader,
        index,
        i,
        storageOffset,
        storageIndex,
        xDims,
        sourceDims
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

int Foam::operator==(const equationOperation& I1, const equationOperation& I2)
{
    return
    (
        I1.sourceType() == I2.sourceType()
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


// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

void Foam::equationOperation::operator=(Foam::equationOperation& eqn)
{
    source_ = eqn.source_;
    sourceIndex_ = eqn.sourceIndex_;
    componentIndex_ = eqn.componentIndex_;
    dictLookupIndex_ = eqn.dictLookupIndex_;
    operation_ = eqn.operation_;
    getSourceScalarFieldFunction_ = eqn.getSourceScalarFieldFunction_;
    opScalarFieldFunction_ = eqn.opScalarFieldFunction_;
    getSourceScalarFunction_ = eqn.getSourceScalarFunction_;
    opScalarFunction_ = eqn.opScalarFunction_;
    getSourceDimsFunction_ = eqn.getSourceDimsFunction_;
    opDimsFunction_ = eqn.opDimsFunction_;
}


// * * * * * * * * * * * * * Friend IOstream Operators * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, equationOperation& I)
{
    label st(I.source_);
    label op(I.operation_);

    is >> st >> I.sourceIndex_ >> I.dictLookupIndex_ >> op;
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const equationOperation& I)
{
    label st(I.source_);
    label op(I.operation_);

    return os   << nl << "/* sourceType:    */\t" << st
                    << " /* "
                    << equationOperation::sourceName(I.source_)
                    << " */" << nl
                << "/* sourceIndex:   */\t" << I.sourceIndex_ << nl
                << "/* componentIndex:*/\t" << I.componentIndex_ << nl
                << "/* dictIndex      */\t" << I. dictLookupIndex_ << nl
                << "/* operation:     */\t" << op
                    << " /* "
                    << equationOperation::opName(I.operation_)
                    << " */" << nl;
}

// ************************************************************************* //
