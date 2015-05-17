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

Foam::dimensionSet Foam::equationReader::internalEvaluateDimensions
(
    const label& equationIndex,
    label storageOffset
) const
{
#   ifdef FULLDEBUG
    // Bounds checking
    if ((equationIndex < 0) || (equationIndex >= size()))
    {
        FatalErrorIn("equationReader::internalEvaluateDimensions")
            << "equationIndex " << equationIndex << " out of bounds (0, "
            << size() - 1 << ")"
            << abort(FatalError);
    }
#   endif

    const equation& eqn(operator[](equationIndex));

    // Launch the reportDimsEvalStartFunction, which does this:
    //  if ((debug == 3) || (debug == 4) || (debug == 5) || (debug == 6))
    //  {
    //      reportDimsEvalStartEnabled(equationIndex);
    //      // reports details to the console that evaluation has commenced
    //  }
    //  else
    //  {
    //      reportDimsEvalStartDisabled(equationIndex);
    //      // does nothing
    //  }
    (*this.*reportDimsEvalStartFunction_)(equationIndex);

    if (eqn.size() == 0)
    {
        parse(equationIndex);
    }

    label storeIndex(-1);

    dimensionSet xDims(dimless);

    for (label i(0); i < eqn.size(); i++)
    {
#       ifdef FULLDEBUG
        if
        (
            (
                (i == 0)
             || (eqn[i - 1].operation() == equationOperation::otstore)
            )
         && (
                eqn[i].operation()
             != equationOperation::otretrieve
            )
        )
        {
            FatalErrorIn("equationReader::internalEvaluateDimensions")
                << "Bad operation list.  Operation at " << i << " either "
                << "follows a 'store', or is the first operation.  Therefore "
                << "it should be retrieve, but it is " << eqn[i].operation()
                << "."
                << abort(FatalError);
        }
#       endif

        // Execute getSource function to which this operation points
        dimensionSet sourceDims
        (
            eqn[i].getSourceDimsFunction
            (
                this,
                equationIndex,
                i,
                storeIndex + storageOffset,
                storageOffset
            )
        );

        // Launch the reportDimsOperationFunction, which does this:
        //  if ((debug == 4) || (debug == 6))
        //  {
        //      reportDimsOperationEnabled(equationIndex, i);
        //      // posts operation-by-operation information to the console
        //  }
        //  else
        //  {
        //      reportDimsOperationDisabled(equationIndex, i);
        //      // does nothing
        //  }
        (*this.*reportDimsOperationFunction_)(equationIndex, i);

        // Execute the eval function to which this operation points
        eqn[i].opDimsFunction
        (
            this,
            equationIndex,
            i,
            storageOffset,
            storeIndex,
            xDims,
            sourceDims
        );

        // Launch the reportDimsResultFunction, which does this:
        //  if ((debug == 4) || (debug == 6))
        //  {
        //      reportDimsResultEnabled(xDims);
        //      // posts result to the console
        //  }
        //  else
        //  {
        //      reportDimsResultDisabled(xDims);
        //      // does nothing
        //  }
        (*this.*reportDimsResultFunction_)(xDims);
    }

    //Move one level back up on the dependents_ list
    if (dependents_.size())
    {
        dependents_.setSize(dependents_.size() - 1);
    }

    storageDims_.setSize(storageOffset);

    // Launch the reportScalarEvalEndFunction, which does this:
    //  if ((debug == 3) || (debug == 4) || (debug == 5) || (debug == 6))
    //  {
    //      reportDimsEvalEndEnabled(xDims);
    //      // reports details to the console that evaluation has completed
    //  }
    //  else
    //  {
    //      reportDimsEvalEndDisabled(xDims);
    //      // does nothing
    //  }
    (*this.*reportDimsEvalEndFunction_)(xDims);

    eqn.setLastResult(xDims);
    return xDims;
}


Foam::scalar Foam::equationReader::internalEvaluateScalar
(
    const label& equationIndex,
    label storageOffset
) const
{
#   ifdef FULLDEBUG
    // Bounds checking
    if ((equationIndex < 0) || (equationIndex >= size()))
    {
        FatalErrorIn("equationReader::internalEvaluateScalar")
            << "equationIndex " << equationIndex << " out of bounds (0, "
            << size() - 1 << ")"
            << abort(FatalError);
    }
#   endif

    const equation& eqn(operator[](equationIndex));

    // Launch the reportScalarEvalStartFunction, which does this:
    //  if ((debug == 1) || (debug == 2) || (debug == 5) || (debug == 6))
    //  {
    //      reportScalarEvalStartEnabled(equationIndex);
    //      // reports details to the console that evaluation has commenced
    //  }
    //  else
    //  {
    //      reportScalarEvalStartDisabled(equationIndex);
    //      // does nothing
    //  }
    (*this.*reportScalarEvalStartFunction_)(equationIndex);

    if (eqn.size() == 0)
    {
        parse(equationIndex);
    }

    label storeIndex(-1);

    scalar x(0.0);

    for (label i(0); i < eqn.size(); i++)
    {
#       ifdef FULLDEBUG
        if
        (
            (
                (i == 0)
             || (eqn[i - 1].operation() == equationOperation::otstore)
            )
         && (
                eqn[i].operation()
             != equationOperation::otretrieve
            )
        )
        {
            FatalErrorIn("equationReader::internalEvaluateScalar")
                << "Bad operation list.  Operation at " << i << " either "
                << "follows a 'store', or is the first operation.  Therefore "
                << "it should be retrieve, but it is " << eqn[i].operation()
                << "."
                << abort(FatalError);
        }
#       endif

        // Execute getSource function to which this operation points
        scalar source
        (
            eqn[i].getSourceScalarFunction
            (
                this,
                equationIndex,
                i,
                storeIndex + storageOffset,
                storageOffset
            )
        );

        // Launch the reportScalarOperationFunction, which does this:
        //  if ((debug == 2) || (debug == 6))
        //  {
        //      reportScalarOperationEnabled(equationIndex, i);
        //      // posts operation-by-operation information to the console
        //  }
        //  else
        //  {
        //      reportScalarOperationDisabled(equationIndex, i);
        //      // does nothing
        //  }
        (*this.*reportScalarOperationFunction_)(equationIndex, i);

        // Execute the eval function to which this operation points
        eqn[i].opScalarFunction
        (
            this,
            equationIndex,
            i,
            storageOffset,
            storeIndex,
            x,
            source
        );

        // Launch the reportScalarResultFunction, which does this:
        //  if ((debug == 2) || (debug == 6))
        //  {
        //      reportScalarResultEnabled(x);
        //      // posts result to the console
        //  }
        //  else
        //  {
        //      reportScalarResultDisabled(x);
        //      // does nothing
        //  }
        (*this.*reportScalarResultFunction_)(x);
    }

    //Move one level back up on the dependents_ list
    if (dependents_.size())
    {
        dependents_.setSize(dependents_.size() - 1);
    }

    storageScalars_.setSize(storageOffset);

    // Launch the reportScalarEvalEndFunction, which does this:
    //  if ((debug == 1) || (debug == 2) || (debug == 5) || (debug == 6))
    //  {
    //      reportScalarEvalEndEnabled(equationIndex);
    //      // reports details to the console that evaluation has completed
    //  }
    //  else
    //  {
    //      reportScalarEvalEndDisabled(equationIndex);
    //      // does nothing
    //  }
    (*this.*reportScalarEvalEndFunction_)(x);

    eqn.setLastResult(x);
    return x;
}


void Foam::equationReader::internalEvaluateScalarField
(
    scalarField& result,
    const label& equationIndex,
    label storageOffset
) const
{
#   ifdef FULLDEBUG
    // Bounds checking
    if ((equationIndex < 0) || (equationIndex >= size()))
    {
        FatalErrorIn("equationReader::internalEvaluateScalarField")
            << "equationIndex " << equationIndex << " out of bounds (0, "
            << size() - 1 << ")"
            << abort(FatalError);
    }
#   endif

    tempSrcField_.setSize(result.size());

    const equation& eqn(operator[](equationIndex));

    // Launch the reportScalarEvalStartFunction, which does this:
    //  if ((debug == 1) || (debug == 2) || (debug == 5) || (debug == 6))
    //  {
    //      reportScalarEvalStartEnabled(equationIndex);
    //      // reports details to the console that evaluation has commenced
    //  }
    //  else
    //  {
    //      reportScalarEvalStartDisabled(equationIndex);
    //      // does nothing
    //  }
    (*this.*reportScalarEvalStartFunction_)(equationIndex);

    if (eqn.size() == 0)
    {
        parse(equationIndex);
    }

    label storeIndex(-1);

    result = 0.0;

    for (label i(0); i < eqn.size(); i++)
    {
#       ifdef FULLDEBUG
        if
        (
            (
                (i == 0)
             || (eqn[i - 1].operation() == equationOperation::otstore)
            )
         && (
                eqn[i].operation()
             != equationOperation::otretrieve
            )
        )
        {
            FatalErrorIn("equationReader::internalEvaluateScalarField")
                << "Bad operation list.  Operation at " << i << " either "
                << "follows a 'store', or is the first operation.  Therefore "
                << "it should be retrieve, but it is " << eqn[i].operation()
                << "."
                << abort(FatalError);
        }
#       endif

        // Execute getSource function to which this operation points
        const scalarField& source
        (
            eqn[i].getSourceScalarFieldFunction
            (
                this,
                equationIndex,
                i,
                storeIndex + storageOffset,
                storageOffset
            )
        );

        // Launch the reportScalarOperationFunction, which does this:
        //  if ((debug == 2) || (debug == 6))
        //  {
        //      reportScalarOperationEnabled(equationIndex, i);
        //      // posts operation-by-operation information to the console
        //  }
        //  else
        //  {
        //      reportScalarOperationDisabled(equationIndex, i);
        //      // does nothing
        //  }
        (*this.*reportScalarOperationFunction_)(equationIndex, i);

        // Execute the eval function to which this operation points
        eqn[i].opScalarFieldFunction
        (
            this,
            equationIndex,
            i,
            storageOffset,
            storeIndex,
            result,
            source
        );

        // Launch the reportScalarResultFunction, which does this:
        //  if ((debug == 2) || (debug == 6))
        //  {
        //      reportScalarResultEnabled(x);
        //      // posts result to the console
        //  }
        //  else
        //  {
        //      reportScalarResultDisabled(x);
        //      // does nothing
        //  }
        (*this.*reportScalarResultFunction_)(result[0]);
    }

    //Move one level back up on the dependents_ list
    if (dependents_.size())
    {
        dependents_.setSize(dependents_.size() - 1);
    }

    storageScalarFields_.setSize(storageOffset);

    // Launch the reportScalarEvalEndFunction, which does this:
    //  if ((debug == 1) || (debug == 2) || (debug == 5) || (debug == 6))
    //  {
    //      reportScalarEvalEndEnabled(equationIndex);
    //      // reports details to the console that evaluation has completed
    //  }
    //  else
    //  {
    //      reportScalarEvalEndDisabled(equationIndex);
    //      // does nothing
    //  }
    (*this.*reportScalarEvalEndFunction_)(result[0]);

    eqn.setLastResult(result[result.size() - 1]);
    tempSrcField_.setSize(0);
}


void Foam::equationReader::checkFinalDimensions
(
    const label& equationIndex,
    dimensionSet& expectedDimensions,
    const word& outputName
) const
{
    const equation& eqn(operator[](equationIndex));
    dimensionSet outputDims(evaluateDimensions(equationIndex));
    if ((outputDims != expectedDimensions) && (dimensionSet::debug))
    {
        WarningIn("equationReader::checkFinalDimenions")
            << "Dimension error thrown for equation " << eqn.name()
            << ".  Output dimensions: " << outputDims << "do not match "
            << "dimensions of destination field " << outputName << ", "
            << expectedDimensions << ". You can initialize " << outputName
            << "'s dimensions with:" << token::NL
            << token::TAB << outputName << ".dimensions().reset" << token::NL
            << token::TAB << "(" << token::NL << token::TAB << token::TAB
            << "equationReaderObject.evaluateDimensions(" << equationIndex
            << "))" << endl;
    }
    expectedDimensions = outputDims;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::equationReader::evaluateScalar
(
    const word& equationName,
    const label cellIndex,
    const label geoIndex
) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateScalar")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    return evaluateScalar(equationIndex, cellIndex, geoIndex);
}


Foam::scalar Foam::equationReader::evaluateScalar
(
    const label equationIndex,
    const label cellIndex,
    const label geoIndex
) const
{
#   ifdef FULLDEBUG
        const equation& eqn(operator[](equationIndex));
        // Index checking
        const labelList& maxFieldSizes(eqn.maxFieldSizes());
        if ((geoIndex < 0) || (cellIndex < 0))
        {
            FatalErrorIn("equationReader::evaluateScalar")
                << "Evaluating " << eqn.name() << ": geoIndex (" << geoIndex
                << ") or cellIndex (" << cellIndex << ") cannot be negative."
                << abort(FatalError);
        }
        else if (maxFieldSizes.size() > 0)
        {
            if (geoIndex >= maxFieldSizes.size())
            {
                FatalErrorIn("equationReader::evaluateScalar")
                    << "Evaluating " << eqn.name() << ": geoIndex ("
                    << geoIndex << ") out of range (0 .. "
                    << maxFieldSizes.size() - 1 << ")."
                    << abort(FatalError);
            }
            else if (cellIndex >= maxFieldSizes[geoIndex])
            {
                FatalErrorIn("equationReader::evaluateScalar")
                    << "Evaluating " << eqn.name() << ": cellIndex ("
                    << cellIndex << ") out of range (0 .. "
                    << maxFieldSizes[geoIndex] - 1 << ")."
                    << abort(FatalError);
            }
        }
#   endif
    cellIndex_ = cellIndex;
    geoIndex_ = geoIndex;
    return internalEvaluateScalar(equationIndex, 0);
}


Foam::dimensionSet
    Foam::equationReader::evaluateDimensions(const word& equationName) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateDimensions")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    return evaluateDimensions(equationIndex);
}


Foam::dimensionSet
    Foam::equationReader::evaluateDimensions(const label equationIndex) const
{
    // Call the associated evaluateDimensions function pointer - has the same
    // effect as this:
    //
    //  const equation& eqn(operator[](equationIndex));
    //  if (eqn.changeDimensions())
    //  {
    //      evaluateDimsDisabled(equationIndex, 0);
    //      // which in turn does: return eqn.overrideDimensions();
    //  }
    //  else
    //  {
    //      evaluadeDimsEnabled(equationIndex);
    //      // which in turn does:
    //      // return internalEvaluateDimensions(equationIndex, 0);
    //  }
    return dimensionSet
    (
        (*this.*evaluateDimsFunctions_[equationIndex])
        (equationIndex, 0)
    );
}


Foam::dimensionedScalar Foam::equationReader::evaluateDimensionedScalar
(
    const word& equationName,
    const label cellIndex,
    const label geoIndex
) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateDimensionedScalar")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    return evaluateDimensionedScalar(equationIndex, cellIndex, geoIndex);
}


Foam::dimensionedScalar Foam::equationReader::evaluateDimensionedScalar
(
    const label equationIndex,
    const label cellIndex,
    const label geoIndex
) const
{
    return dimensionedScalar
    (
        operator[](equationIndex).name(),
        evaluateDimensions(equationIndex),
        evaluateScalar(equationIndex, cellIndex, geoIndex)
    );
}


void Foam::equationReader::evaluateScalarField
(
    scalarField& resultField,
    const word& equationName,
    const label geoIndex
) const
{
    label equationIndex(lookup(equationName));
    if (equationIndex < 0)
    {
        FatalErrorIn("equationReader::evaluateScalarField")
            << "Equation name " << equationName << " not found."
            << abort(FatalError);
    }
    evaluateScalarField(resultField, equationIndex, geoIndex);
}


void Foam::equationReader::evaluateScalarField
(
    scalarField& resultField,
    const label equationIndex,
    const label geoIndex
) const
{
#   ifdef FULLDEBUG
        const equation& eqn(operator[](equationIndex));
        // Index checking
        const labelList& maxFieldSizes(eqn.maxFieldSizes());
        if (geoIndex < 0)
        {
            FatalErrorIn("equationReader::evaluateScalarField")
                << "Evaluating " << eqn.name() << ": geoIndex (" << geoIndex
                << ") cannot be negative."
                << abort(FatalError);
        }
        else if (maxFieldSizes.size() > 0)
        {
            if (geoIndex >= maxFieldSizes.size())
            {
                FatalErrorIn("equationReader::evaluateScalarField")
                    << "Evaluating " << eqn.name() << ": geoIndex ("
                    << geoIndex << ") out of range (0 .. "
                    << maxFieldSizes.size() - 1 << ")."
                    << abort(FatalError);
            }
            else if (resultField.size() != maxFieldSizes[geoIndex])
            {
                FatalErrorIn("equationReader::evaluateScalarField")
                    << "Evaluating " << eqn.name() << ": field size mismatch. "
                    << "result.size() = " << resultField.size() << ", "
                    << "expected = " << maxFieldSizes[geoIndex] - 1 << "."
                    << abort(FatalError);
            }
        }
#   endif
    geoIndex_ = geoIndex;
    if (resultField.size())
    {
        internalEvaluateScalarField(resultField, equationIndex, 0);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
