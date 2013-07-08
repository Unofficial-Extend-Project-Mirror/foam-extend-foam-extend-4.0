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

Foam::dimensionedScalar Foam::equationReader::getSourceNone
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    return dimensionedScalar("noSource", dimless, 0);
}


Foam::dimensionedScalar Foam::equationReader::getSourceDictSourceDScalar
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    word varName(dictLookups_[eqOp.dictLookupIndex()]);
    
    ITstream srcStrm
    (
        dictSources_[zeroSourceIndex].lookup(varName)
    );
    srcStrm >> returnMe;
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceDictSourceScalar
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    word varName(dictLookups_[eqOp.dictLookupIndex()]);
    
    returnMe.name() = varName;
    returnMe.value() = readScalar
    (
        dictSources_[zeroSourceIndex].lookup(varName)
    );
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceExternalDScalar
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    dsEqual(returnMe, externalDScalars_[zeroSourceIndex]);
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceExternalScalar
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    returnMe.name() = externalScalarNames_[zeroSourceIndex];
    returnMe.value() = externalScalars_[zeroSourceIndex];
    returnMe.dimensions().reset
    (
        externalScalarDimensions_[zeroSourceIndex]
    );
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceExternalScalarList
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    dsEqual
    (
        returnMe,
        dimensionedScalar
        (
            externalScalarListNames_[zeroSourceIndex],
            externalScalarListDimensions_[zeroSourceIndex],
            externalScalarLists_
                [zeroSourceIndex]
                [externalScalarListIndex_[zeroSourceIndex]]
        )
    );
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceInternalScalar
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    returnMe.name() = "internalConstant";
    returnMe.value() = internalScalars_[zeroSourceIndex];
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceEquation
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
    const equationOperation& eqOp(eqn[equationOperationIndex]);
    label zeroSourceIndex = mag(eqOp.sourceIndex()) - 1;

    dependents_.setSize(dependents_.size() + 1);
    dependents_[dependents_.size() - 1] = equationIndex;
    if (debug)
    {
        Info << "Embedded equation dispatch." << endl;
    }
    dsEqual(returnMe, evaluate(zeroSourceIndex, maxStoreIndex + 1));
    if (debug)
    {
        Info << "Returned from embedded equation." << endl;
    }
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceEquationCircRefDetect
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    equation& eqn(eqns_[equationIndex]);
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
                    eqns_[j].equationName()
                );
                dependencies.append("-->");
            }
            dependencies.append(eqns_[i].equationName());
            FatalErrorIn("equationReader::getSource")
                << "Circular reference detected when evaluating "
                << "the equation for " << eqn.equationName()
                << ", given by:" << token::NL << token::TAB
                << eqn.rawText() << token::NL << "The circular "
                << "dependency is:" << token::NL << token::TAB
                << dependencies
                << abort(FatalError);
        }
    }
    if (debug)
    {
        Info << "Embedded equation dispatch." << endl;
    }
    dsEqual(returnMe, evaluate(zeroSourceIndex, maxStoreIndex + 1));
    eqOp.assignSourceFunction
    (
        &Foam::equationReader::getSourceEquation
    );
    if (debug)
    {
        Info << "Returned from embedded equation." << endl;
    }
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


Foam::dimensionedScalar Foam::equationReader::getSourceStorage
(
    equationReader * eqnReader,
    const label equationIndex,
    const label equationOperationIndex,
    const label maxStoreIndex,
    const label storageOffset
)
{
    dimensionedScalar returnMe("noSource", dimless, 0);
    const equation& eqn(eqns_[equationIndex]);
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
    dsEqual(returnMe, storage_[zeroSourceIndex + storageOffset]);
    returnMe.value() = sign(eqOp.sourceIndex()) * returnMe.value();
    return returnMe;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
