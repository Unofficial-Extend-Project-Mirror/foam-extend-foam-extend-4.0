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

Foam::scalar Foam::equationReader::getScalarSrcNone
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    return 0.0;
}


Foam::scalar Foam::equationReader::getScalarSrcStorage
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
        FatalErrorIn("equationReader::getSouce")
            << "Index " << zeroSourceIndex << " out of bounds (0, "
            << maxStoreIndex - storageOffset << ")"
            << abort(FatalError);
    }
#   endif
    scalar returnMe(storageScalars_[zeroSourceIndex + storageOffset]);
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcActiveSource
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

    scalar returnMe
    (
        activeSources_[zeroSourceIndex].evaluateScalar
        (
            eqOp.componentIndex(),
            cellIndex_,
            geoIndex_
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcEquation
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

    scalar returnMe
    (
        internalEvaluateScalar(zeroSourceIndex, maxStoreIndex + 1)
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

    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcEquationCircRefDetect
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
            FatalErrorIn("equationReader::getScalarSrcEquationCircRefDetect")
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
    scalar returnMe
    (
        internalEvaluateScalar(zeroSourceIndex, maxStoreIndex + 1)
    );
    eqOp.assignSourceScalarFunction
    (
        &Foam::equationReader::getScalarSrcEquation
    );
    eqOp.assignSourceScalarFieldFunction
    (
        &Foam::equationReader::getScalarFieldSrcEquation
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

    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcInternalScalar
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

    scalar returnMe(internalScalars_[zeroSourceIndex]);
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcDictSourceDScalar
(
    const equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
) const
{
    dimensionedScalar ds("noSource", dimless, 0.0);
    const equation& eqn(operator[](equationIndex));
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    word varName(dictLookups_[eqOp.dictLookupIndex()]);

    ITstream srcStrm
    (
        dictSources_[zeroSourceIndex].lookup(varName)
    );
    srcStrm >> ds;
    return ds.value() * sign(eqOp.sourceIndex());
}


Foam::scalar Foam::equationReader::getScalarSrcDictSourceScalar
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

    word varName(dictLookups_[eqOp.dictLookupIndex()]);

    scalar returnMe
    (
        readScalar(dictSources_[zeroSourceIndex].lookup(varName))
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcScalarSource
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
    scalar returnMe
    (
        scalarSources_.singleValue
        (
            zeroSourceIndex,
            eqOp.componentIndex()
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcScalarFieldSource
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
    scalar returnMe
    (
        scalarSources_.fieldValue
        (
            zeroSourceIndex,
            eqOp.componentIndex(),
            cellIndex_,
            geoIndex_
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcVectorSource
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
    scalar returnMe
    (
        vectorSources_.singleValue
        (
            zeroSourceIndex,
            eqOp.componentIndex()
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcVectorFieldSource
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
    scalar returnMe
    (
        vectorSources_.fieldValue
        (
            zeroSourceIndex,
            eqOp.componentIndex(),
            cellIndex_,
            geoIndex_
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcTensorSource
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
    scalar returnMe
    (
        tensorSources_.singleValue
        (
            zeroSourceIndex,
            eqOp.componentIndex()
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcTensorFieldSource
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
    scalar returnMe
    (
        tensorSources_.fieldValue
        (
            zeroSourceIndex,
            eqOp.componentIndex(),
            cellIndex_,
            geoIndex_
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcDiagTensorSource
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
    scalar returnMe
    (
        diagTensorSources_.singleValue
        (
            zeroSourceIndex,
            eqOp.componentIndex()
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcDiagTensorFieldSource
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
    scalar returnMe
    (
        diagTensorSources_.fieldValue
        (
            zeroSourceIndex,
            eqOp.componentIndex(),
            cellIndex_,
            geoIndex_
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcSymmTensorSource
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
    scalar returnMe
    (
        symmTensorSources_.singleValue
        (
            zeroSourceIndex,
            eqOp.componentIndex()
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcSymmTensorFieldSource
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
    scalar returnMe
    (
        symmTensorSources_.fieldValue
        (
            zeroSourceIndex,
            eqOp.componentIndex(),
            cellIndex_,
            geoIndex_
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcSphericalTensorSource
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
    scalar returnMe
    (
        sphericalTensorSources_.singleValue
        (
            zeroSourceIndex,
            eqOp.componentIndex()
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}


Foam::scalar Foam::equationReader::getScalarSrcSphericalTensorFieldSource
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
    scalar returnMe
    (
        sphericalTensorSources_.fieldValue
        (
            zeroSourceIndex,
            eqOp.componentIndex(),
            cellIndex_,
            geoIndex_
        )
    );
    returnMe *= sign(eqOp.sourceIndex());
    return returnMe;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
