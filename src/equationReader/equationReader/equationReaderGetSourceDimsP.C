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

Foam::dimensionSet Foam::equationReader::getDimsSrcNone
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    return dimless;
}


Foam::dimensionSet Foam::equationReader::getDimsSrcStorage
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

#   ifdef FULLDEBUG
    if ((zeroSourceIndex + storageOffset) > maxStoreIndex)
    {
        FatalErrorIn("equationReader::getDimsSrcStorage")
            << "Index " << zeroSourceIndex << " out of bounds (0, "
            << maxStoreIndex - storageOffset << ")"
            << abort(FatalError);
    }
#   endif
    return storageDims_[zeroSourceIndex + storageOffset];
}


Foam::dimensionSet Foam::equationReader::getDimsSrcActiveSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return activeSources_[zeroSourceIndex].dimensions();
}


Foam::dimensionSet Foam::equationReader::getDimsSrcEquation
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    dependents_.setSize(dependents_.size() + 1);
    dependents_[dependents_.size() - 1] = equationIndex;

    // Launch the reportEmbeddedDispatchFunction:
    //  if (debug)
    //  {
    //      reportEmbeddedDispatchEnabled;
    //      // or: Info << "Embedded equation dispatch." << endl;
    //  }
    //  else
    //  {
    //      reportEmbeddedDispatchDisabled();
    //      // does nothing
    //  }
    (*this.*reportEmbeddedDispatchFunction_)();

    // Call the associated evaluateDimensions function pointer - has the same
    // effect as this:
    //
    //  const equation& eqn(operator[](zeroSourceIndex));
    //  if (eqn.changeDimensions())
    //  {
    //      evaluateDimsDisabled(zeroSourceIndex, maxStoreIndex + 1);
    //      // which in turn does: return eqn.overrideDimensions();
    //  }
    //  else
    //  {
    //      evaluadeDimsEnabled(zeroSourceIndex, maxStoreIndex + 1);
    //      // which in turn does:
    //      // return internalEvaluateDimensions
    //      // (zeroSourceIndex, maxStoreIndex + 1);
    //  }
    dimensionSet returnMe
    (
        (*this.*evaluateDimsFunctions_[zeroSourceIndex])
        (zeroSourceIndex, maxStoreIndex + 1)
    );

    // Launch the reportEmbeddedReturnFunction:
    //  if (debug)
    //  {
    //      reportEmbeddedReturnEnabled;
    //      // or: Info << "Return from equation equation." << endl;
    //  }
    //  else
    //  {
    //      reportEmbeddedReturnDisabled();
    //      // does nothing
    //  }
    (*this.*reportEmbeddedReturnFunction_)();

    //Move one level back up on the dependents_ list
    if (dependents_.size())
    {
        dependents_.setSize(dependents_.size() - 1);
    }

    return returnMe;
}


Foam::dimensionSet Foam::equationReader::getDimsSrcEquationCircRefDetect
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    // Check for circular references
    dependents_.setSize(dependents_.size() + 1);
    dependents_[dependents_.size() - 1] = equationIndex;
    forAll(dependents_, i)
    {
        if (dependents_[i] == zeroSourceIndex)
        {
            // Circular reference detected
            string dependencies;
            for (label j(i); j < dependents_.size(); j++)
            {
                dependencies.append
                (
                    operator[](dependents_[j]).name()
                );
                dependencies.append("-->");
            }
            dependencies.append(operator[](dependents_[i]).name());
            FatalErrorIn
            (
                "equationReader::getDimsSrcEquationCircRefDetect"
            )
                << "Circular reference detected when evaluating "
                << "the equation for " << eqn.name()
                << ", given by:" << token::NL << token::TAB
                << eqn.rawText() << token::NL << "The circular "
                << "dependency is:" << token::NL << token::TAB
                << dependencies
                << abort(FatalError);
        }
    }

    // Launch the reportEmbeddedDispatchFunction:
    //  if (debug)
    //  {
    //      reportEmbeddedDispatchEnabled;
    //      // or: Info << "Embedded equation dispatch." << endl;
    //  }
    //  else
    //  {
    //      reportEmbeddedDispatchDisabled();
    //      // does nothing
    //  }
    (*this.*reportEmbeddedDispatchFunction_)();

    // Call the associated evaluateDimensions function pointer - has the same
    // effect as this:
    //
    //  const equation& eqn(operator[](zeroSourceIndex));
    //  if (eqn.changeDimensions())
    //  {
    //      evaluateDimsDisabled(zeroSourceIndex, maxStoreIndex + 1);
    //      // which in turn does: return eqn.overrideDimensions();
    //  }
    //  else
    //  {
    //      evaluadeDimsEnabled(zeroSourceIndex, maxStoreIndex + 1);
    //      // which in turn does:
    //      // return internalEvaluateDimensions
    //      // (zeroSourceIndex, maxStoreIndex + 1);
    //  }
    dimensionSet returnMe
    (
        (*this.*evaluateDimsFunctions_[zeroSourceIndex])
        (zeroSourceIndex, maxStoreIndex + 1)
    );

    // Launch the reportEmbeddedReturnFunction:
    //  if (debug)
    //  {
    //      reportEmbeddedReturnEnabled;
    //      // or: Info << "Return from equation equation." << endl;
    //  }
    //  else
    //  {
    //      reportEmbeddedReturnDisabled();
    //      // does nothing
    //  }
    (*this.*reportEmbeddedReturnFunction_)();

    //Move one level back up on the dependents_ list
    if (dependents_.size())
    {
        dependents_.setSize(dependents_.size() - 1);
    }

    return returnMe;
}


Foam::dimensionSet Foam::equationReader::getDimsSrcInternalScalar
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    return dimless;
}


Foam::dimensionSet Foam::equationReader::getDimsSrcDictSourceDScalar
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    dimensionedScalar ds("noSource", dimless, 0);
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    word varName(dictLookups_[eqOp.dictLookupIndex()]);

    ITstream srcStrm
    (
        dictSources_[zeroSourceIndex].lookup(varName)
    );
    srcStrm >> ds;
    return ds.dimensions();
}


Foam::dimensionSet Foam::equationReader::getDimsSrcDictSourceScalar
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    return dimless;
}


Foam::dimensionSet Foam::equationReader::getDimsSrcScalarSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return scalarSources_.singleDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcScalarFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return scalarSources_.fieldDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcVectorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return vectorSources_.singleDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcVectorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return vectorSources_.fieldDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return tensorSources_.singleDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return tensorSources_.fieldDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcDiagTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return diagTensorSources_.singleDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcDiagTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return diagTensorSources_.fieldDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcSymmTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return symmTensorSources_.singleDimensions(zeroSourceIndex);
}

Foam::dimensionSet Foam::equationReader::getDimsSrcSymmTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return symmTensorSources_.fieldDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcSphericalTensorSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return sphericalTensorSources_.singleDimensions(zeroSourceIndex);
}


Foam::dimensionSet Foam::equationReader::getDimsSrcSphericalTensorFieldSource
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;
    return sphericalTensorSources_.fieldDimensions(zeroSourceIndex);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
