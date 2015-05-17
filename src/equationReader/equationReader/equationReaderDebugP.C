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

void Foam::equationReader::reportEmbeddedDispatchDisabled () const
{
    // do nothing
}


void Foam::equationReader::reportEmbeddedDispatchEnabled () const
{
    Info << "Embedded equation dispatch." << endl;
}


void Foam::equationReader::reportEmbeddedReturnDisabled () const
{
    // do nothing
}


void Foam::equationReader::reportEmbeddedReturnEnabled () const
{
    Info << "Returned from embedded equation." << endl;
}


void Foam::equationReader::reportScalarEvalStartDisabled
(
    const label& index
) const
{
    // Do nothing
}


void Foam::equationReader::reportScalarEvalStartEnabled
(
    const label& equationIndex
) const
{
    const equation& eqn(operator[](equationIndex));
    Info << "Evaluating equation " << equationIndex << ", "
        << eqn.name() << " at (geoIndex, cellIndex)=("
        << geoIndex_ << ", " << cellIndex_ << "), given by:"
        << token::NL << token::TAB << eqn.rawText() << endl;
}


void Foam::equationReader::reportScalarOperationDisabled
(
    const label& index,
    const label& i
) const
{
    // do nothing
}


void Foam::equationReader::reportScalarOperationEnabled
(
    const label& index,
    const label& i
) const
{
    const equationOperation& eqop(operator[](index)[i]);
    const label zeroSourceIndex(mag(eqop.sourceIndex()) - 1);
    Info << "Performing operation: ["
        << equationOperation::opName(eqop.operation()) << "] using source [";
    switch (eqop.sourceType())
    {
        case equationOperation::stnone:
            Info << "none";
            break;
        case equationOperation::ststorage:
            Info << "memory spot (" << zeroSourceIndex << ")";
            break;
        case equationOperation::stactiveSource:
            Info << activeSourceNames_[zeroSourceIndex];
            break;
        case equationOperation::stequation:
            Info << operator[](zeroSourceIndex).name();
            break;
        case equationOperation::stinternalScalar:
            Info << "constant (" << internalScalars_[zeroSourceIndex] << ")";
            break;
        case equationOperation::stdictSource:
            Info << dictLookups_[zeroSourceIndex];
            break;
        case equationOperation::stscalarSource:
            Info << scalarSources_.singleName(zeroSourceIndex);
            break;
        case equationOperation::stscalarFieldSource:
            Info << scalarSources_.fieldName(zeroSourceIndex);
            break;
        case equationOperation::stvectorSource:
            Info << vectorSources_.singleName(zeroSourceIndex);
            break;
        case equationOperation::stvectorFieldSource:
            Info << vectorSources_.fieldName(zeroSourceIndex);
            break;
        case equationOperation::sttensorSource:
            Info << tensorSources_.singleName(zeroSourceIndex);
            break;
        case equationOperation::sttensorFieldSource:
            Info << tensorSources_.fieldName(zeroSourceIndex);
            break;
        case equationOperation::stdiagTensorSource:
            Info << diagTensorSources_.singleName(zeroSourceIndex);
            break;
        case equationOperation::stdiagTensorFieldSource:
            Info << diagTensorSources_.fieldName(zeroSourceIndex);
            break;
        case equationOperation::stsymmTensorSource:
            Info << symmTensorSources_.singleName(zeroSourceIndex);
            break;
        case equationOperation::stsymmTensorFieldSource:
            Info << symmTensorSources_.fieldName(zeroSourceIndex);
            break;
        case equationOperation::stsphericalTensorSource:
            Info << sphericalTensorSources_.singleName(zeroSourceIndex);
            break;
        case equationOperation::stsphericalTensorFieldSource:
            Info << sphericalTensorSources_.fieldName(zeroSourceIndex);
            break;
    }
    Info << "] read from ["
        << equationOperation::sourceName(eqop.sourceType()) << "]..." << endl;
}


void Foam::equationReader::reportScalarResultDisabled
(
    const scalar& x
) const
{
    // do nothing
}


void Foam::equationReader::reportScalarResultEnabled
(
    const scalar& x
) const
{
    Info << "Operation result is " << x << endl;
}


void Foam::equationReader::reportScalarEvalEndDisabled
(
    const scalar& x
) const
{
    // Do nothing
}


void Foam::equationReader::reportScalarEvalEndEnabled
(
    const scalar& x
) const
{
    Info << "Equation evaluated.  Result is: " << x << endl;
}


void Foam::equationReader::reportDimsEvalStartDisabled
(
    const label& index
) const
{
    // Do nothing
}


void Foam::equationReader::reportDimsEvalStartEnabled
(
    const label& equationIndex
) const
{
    const equation& eqn(operator[](equationIndex));
    Info << "Evaluating equation " << equationIndex << ", "
        << eqn.name() << " at (geoIndex, cellIndex)=("
        << geoIndex_ << ", " << cellIndex_ << "), given by:"
        << token::NL << token::TAB << eqn.rawText() << endl;
}


void Foam::equationReader::reportDimsOperationDisabled
(
    const label& index,
    const label& i
) const
{
    // do nothing
}


void Foam::equationReader::reportDimsOperationEnabled
(
    const label& index,
    const label& i
) const
{
    reportScalarOperationEnabled(index, i);
}


void Foam::equationReader::reportDimsResultDisabled
(
    const dimensionSet& xDims
) const
{
    // do nothing
}


void Foam::equationReader::reportDimsResultEnabled
(
    const dimensionSet& xDims
) const
{
    Info << "Operation result is " << xDims << endl;
}


void Foam::equationReader::reportDimsEvalEndDisabled
(
    const dimensionSet& xDims
) const
{
    // Do nothing
}


void Foam::equationReader::reportDimsEvalEndEnabled
(
    const dimensionSet& xDims
) const
{
    Info << "Equation evaluated.  Result is: " << xDims << endl;
}

// ************************************************************************* //
