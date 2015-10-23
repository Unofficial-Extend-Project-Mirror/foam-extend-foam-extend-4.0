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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::labelList Foam::multiSolver::findSuperLoops(const fileName& path)
{
    fileNameList dirEntries(readDir(path, fileName::DIRECTORY));

    labelList superLoopList(dirEntries.size());
    label nSuperLoops(0);

    // Loop through dirEntries, checking for valid integers, sort entries
    forAll(dirEntries, de)
    {
        // Check if directory is "initial"
        if (dirEntries[de] == "initial")
        {
            superLoopList[nSuperLoops++] = -1;
            continue;
        }

        IStringStream superLoopStream(dirEntries[de]);
        token superLoopToken(superLoopStream);

        // Check if directory is an integer
        if (superLoopToken.isLabel() && superLoopStream.eof())
        {
            superLoopList[nSuperLoops++] = superLoopToken.labelToken();
        }
    }
    superLoopList.setSize(nSuperLoops);
    sort(superLoopList);
    return superLoopList;
}


Foam::timeCluster Foam::multiSolver::findClosestGlobalTime
(
    const Foam::scalar value,
    const Foam::timeClusterList& tcl,
    const bool& exact
)
{
    scalar maxDiff(-VGREAT);
    scalar minDiff(VGREAT);
    label underShoot(-1);
    label best(-1);
    label initial(-1);

    // Find closest global minimum that does not exceed value
    forAll(tcl, i)
    {
        if (!tcl[i].times().size()) continue;
        scalar diff(tcl[i].globalMinValue() - value);
        if ((diff <= 0) && (diff >= maxDiff))
        {
            // "initial" directory may be a duplicate match - others take
            // priority.
            if ((tcl[i].superLoop() < 0) && (diff > maxDiff))
            {
                initial = i;
            }
            else
            {
                initial = -1;
                best = i;
            }
            maxDiff = diff;
        }
        else if ((diff > 0) && (diff < minDiff))
        {
            // This is in case all timeClusters exceed value - then we have to
            // return the closest minValue
            minDiff = diff;
            underShoot = i;
        }
    }

    if (initial != -1)
    {
        // "initial" directory is the only match
        best = initial;
    }

    if (best == -1)
    {
        if (minDiff != -1)
        {
            // All timeClusters exceed value, return closest minValue
            best = underShoot;
        }
        else
        {

            FatalErrorIn("multiSolver::findClosestGlobalTime")
                << "The timeClusterList passed to this function has no non-"
                << "empty instantLists.  Use timeClusterList::purgeEmpties "
                << "and check its return value to prevent this."
                << abort(FatalError);
        }
    }

    label timeIndex
    (
        Time::findClosestTimeIndex
        (
            tcl[best].times(), value - tcl[best].globalOffset()
        )
    );

    if (exact && (maxDiff < -VSMALL))
    {
        FatalErrorIn("multiSolver::findClosestGlobalTime")
            << "No exact match found for global time = " << value
            << abort(FatalError);
    }

    return tcl[best](timeIndex);
}


Foam::timeCluster Foam::multiSolver::findClosestLocalTime
(
    const Foam::scalar value,
    const Foam::timeClusterList& tcl,
    const bool& exact
)
{
    scalar maxDiff(-VGREAT);
    scalar minDiff(VGREAT);
    label underShoot(-1);
    label best(-1);
    label initial(-1);
    timeClusterList tclDummy;
    const timeClusterList * tclPtr;

    if (nonOverlapping(tcl))
    {
        tclPtr = & tcl;
    }
    else
    {
        tclDummy = tcl.selectiveSubList(findMaxSuperLoopIndices(tcl));
        tclPtr = & tclDummy;
    }

    for (label i = 0; i < tclPtr->size(); i++)
    {
        if (!tclPtr->operator[](i).times().size()) continue;
        scalar diff(tclPtr->operator[](i).localMinValue() - value);

        if ((diff <= 0) && (diff >= maxDiff))
        {
            // "initial" directory may be a duplicate match - others take
            // priority.

            if ((tclPtr->operator[](i).superLoop() < 0) && (diff > maxDiff))
            {
                initial = i;
            }
            else
            {
                initial = -1;
                best = i;
            }
            maxDiff = diff;
        }
        else if ((diff > 0) && (diff < minDiff))
        {
            // This is in case all timeClusters exceed value - then we have to
            // return the closest minValue
            minDiff = diff;
            underShoot = i;
        }
    }

    if (initial != -1)
    {
        // "initial" directory is the only match
        best = initial;
    }

    if (best == -1)
    {
        if (minDiff != -1)
        {
            // All timeClusters exceed value, return closest minValue
            best = underShoot;
        }
        else
        {

            FatalErrorIn("multiSolver::findClosestLocalTime")
                << "The timeClusterList passed to this function has no non-"
                << "empty instantLists.  Use timeClusterList::purgeEmpties "
                << "and check its return value to prevent this."
                << abort(FatalError);
        }
    }

    label timeIndex
    (
        Time::findClosestTimeIndex
        (
            tclPtr->operator[](best).times(), value
        )
    );

    if (exact && (maxDiff < -VSMALL))
    {
        FatalErrorIn("multiSolver::findClosestLocalTime")
            << "No exact match found for local time = " << value
            << abort(FatalError);
    }

    return tclPtr->operator[](best)(timeIndex);
}


Foam::timeCluster Foam::multiSolver::findLatestGlobalTime
(
    const Foam::timeClusterList& tcl
)
{
    timeCluster bestMax(0);

    timeCluster currentMax;

    forAll(tcl, i)
    {
        if (tcl[i].times().size() == 0) continue;

        currentMax = tcl[i](tcl[i].globalMaxIndex());
        if
        (
            (currentMax.globalValue(0) > bestMax.globalValue(0))
         || (
                (currentMax.globalValue(0) == bestMax.globalValue(0))
             && (currentMax.superLoop() != -1)
            )
        )
        {
            bestMax = currentMax;
        }
    }

    if (bestMax.solverDomainName() == word::null)
    {
        FatalErrorIn("multiSolver::findLatestGlobalTime")
            << "The timeClusterList passed to this function has no non-empty "
            << "instantLists.  Use timeClusterList::purgeEmpties and check its"
            << " return value to prevent this."
            << abort(FatalError);
    }
    return bestMax;
}


Foam::timeCluster Foam::multiSolver::findLatestLocalTime
(
    const Foam::timeClusterList& tcl
)
{
    timeClusterList dummyTcl;
    const timeClusterList * tclPtr;
    timeCluster bestMax(0);
    timeCluster currentMax;

    if (nonOverlapping(tcl))
    {
        tclPtr = & tcl;
    }
    else
    {
        dummyTcl = tcl.selectiveSubList(findMaxSuperLoopIndices(tcl));
        tclPtr = & dummyTcl;
    }

    for (label i = 0; i < tclPtr->size(); i++)
    {
        if (tclPtr->operator[](i).times().size() == 0) continue;

        currentMax =
            tclPtr->operator[](i)(tclPtr->operator[](i).localMaxIndex());

        if
        (
            (currentMax.localValue(0) > bestMax.localValue(0))
         || (
                (currentMax.localValue(0) == bestMax.localValue(0))
             && (currentMax.superLoop() != -1)
            )
        )
        {
            bestMax = currentMax;
        }
    }

    if (bestMax.solverDomainName() == word::null)
    {
        FatalErrorIn("multiSolver::findLatestLocalTime")
            << "The timeClusterList passed to this function has no non-empty "
            << "instantLists.  Use timeClusterList::purgeEmpties and check its"
            << " return value to prevent this."
            << abort(FatalError);
    }

    return bestMax;
}


Foam::timeCluster Foam::multiSolver::findGlobalIndex
(
    const label& index,
    const timeClusterList& tcl
)
{
    forAll(tcl, i)
    {
        if (tcl[i].globalIndex() == index)
        {
            return tcl[i];
        }
    }
    return timeCluster();
}


void Foam::multiSolver::includePreviousTimes
(
    timeCluster& tc
) const
{
    scalar minTime(VGREAT);
    forAll(tc, i)
    {
        minTime = min(tc[i].value(), minTime);
    }
    if (minTime != VGREAT)
    {
        fileName thePath(findInstancePath(tc, 0).path());
        instantList il(Time::findTimes(thePath));
        forAll(il, i)
        {
            if (il[i].value() < minTime)
            {
                label newIndex(tc.size());
                tc.setSize(newIndex + 1);
                tc[newIndex] = il[i];
            }
        }
        if (tc.size() > 1)
        {
            std::sort(&tc[0], tc.end(), instant::less());
        }
    }
}


Foam::fileName Foam::multiSolver::findInstancePath
(
    const timeCluster& tc,
    const label& index
) const
{
    if (!tc.times().size())
    {
        FatalErrorIn("multiSolver::findInstancePath")
            << "The timeClusterList passed to this function has no non-empty "
            << "instantLists.  Use timeClusterList::purgeEmpties and check its"
            << " return value to prevent this."
            << abort(FatalError);
    }
    if (tc.superLoop() < 0)
    {
        // Initial directory
        return fileName
        (
            multiDictRegistry_.path()/"multiSolver"/tc.solverDomainName()
                /"initial"/tc[index].name()
        );
    }
    else
    {
        return fileName
        (
            multiDictRegistry_.path()/"multiSolver"/tc.solverDomainName()
                /name(tc.superLoop())/tc[index].name()
        );
    }
}


Foam::label
    Foam::multiSolver::findMaxSuperLoopValue(const timeClusterList& tcl)
{
    if (!tcl.size())
    {
        FatalErrorIn("multiSolver::findMaxSuperLoopValue")
            << "The timeClusterList passed to this function is empty.  Use "
            << "timeClusterList::purgeEmpties and check its return value to "
            << "prevent this."
            << abort(FatalError);
    }
    return tcl[findMaxSuperLoopIndices(tcl)[0]].superLoop();
}


Foam::labelList
    Foam::multiSolver::findMaxSuperLoopIndices(const timeClusterList& tcl)
{
    label currentMax(-2);
    labelList bestIndices(0);
    label maxesFound(0);

    forAll(tcl, i)
    {
        if (!tcl[i].times().size()) continue;
        if (currentMax == tcl[i].superLoop())
        {
            bestIndices.setSize(++maxesFound);
            bestIndices[maxesFound - 1] = i;
        }
        else if (currentMax < tcl[i].superLoop())
        {
            currentMax = tcl[i].superLoop();
            maxesFound = 1;
            bestIndices.setSize(1);
            bestIndices[0] = i;
        }
    }

    if (bestIndices.size() == 0)
    {
        FatalErrorIn("multiSolver::findMaxSuperLoopIndices")
            << "The timeClusterList passed to this function is empty.  Use "
            << "timeClusterList::purgeEmpties and check its return value to "
            << "prevent this."
            << abort(FatalError);
    }
    return bestIndices;
}


bool Foam::multiSolver::nonOverlapping
(

    const Foam::timeClusterList& tcl,
    const bool useGlobalTime
)
{
    // see tuple2Lists.H
    scalarScalarList range(tcl.size());

    if (useGlobalTime)
    {
        forAll(tcl, i)
        {
            if (!tcl[i].times().size())
            {
                range[i].first() = 0;
                range[i].second() = 0;
            }
            else
            {
                range[i].first() = tcl[i].globalMinValue();
                range[i].second() = tcl[i].globalMaxValue();
            }
        }
    }
    else
    {
        forAll(tcl, i)
        {
            if (!tcl[i].times().size())
            {
                range[i].first() = 0;
                range[i].second() = 0;
            }
            else
            {
                range[i].first() = tcl[i].localMinValue();
                range[i].second() = tcl[i].localMaxValue();
            }
        }
    }
    sortTuple2ListBy1stThen2nd(range);

    for (label i = 0; i < (range.size() - 1); i++)
    {
        // Using '-SMALL' below is a temporary bug fix
        if (range[i + 1].first() - range[i].second() < -SMALL)
        {
            // timeClusters overlap
            return false;
        }
    }
    return true;
}


Foam::timeCluster Foam::multiSolver::readSuperLoopTimes
(
    const Foam::word& solverDomain,
    const Foam::label superLoop,
    const Foam::word& processor
) const
{
    fileName currentPath;
    if (processor.size())
    {
        currentPath = multiDictRegistry_.path()/processor;
    }
    else
    {
        currentPath = multiDictRegistry_.path();
    }
    if (superLoop < 0)
    {
        currentPath = currentPath/"multiSolver"/solverDomain
            /"initial";
    }
    else
    {
        currentPath = currentPath/"multiSolver"/solverDomain
            /name(superLoop);
    }

    fileName mstFileName
    (
        currentPath/"multiSolverTime"
    );
    IFstream mstFile(mstFileName);

    bool mstFileGood(false);
    scalar globalOffset(0);
    label globalIndex(0);

    if (mstFile.good())
    {
        dictionary mstDict(mstFile);
        if (mstDict.found("globalOffset"))
        {
            globalOffset =
                readScalar(mstDict.lookup("globalOffset"));
            globalIndex =
                readLabel(mstDict.lookup("globalIndex"));
            mstFileGood = true;
        }
    }

    if ((!mstFileGood) && (superLoop != -1))
    {
        WarningIn("multiSolver::readSuperLoopTimes")
            << "Bad or missing multiSolverTime dictionary (auto-"
            << "generated) in case/multiSolver/" << solverDomain
            << "/" << superLoop << ".  Assuming globalOffset = 0"
            << endl;
    }
    timeCluster tc
    (
        Time::findTimes(currentPath),
        globalOffset,
        globalIndex,
        superLoop,
        solverDomain
    );
    return tc;
}


Foam::timeClusterList Foam::multiSolver::readSolverDomainTimes
(
    const word& solverDomain,
    const word processor
) const
{
    timeClusterList tcl(0);
    label nTimeClusters(0);

    fileName locale;
    if (processor.size())
    {
        locale = processor/"multiSolver";
    }
    else
    {
        locale = "multiSolver";
    }
    fileName currentPath
    (
        multiDictRegistry_.path()/locale/solverDomain
    );

    labelList superLoopList(multiSolver::findSuperLoops(currentPath));

    // Loop through superLoopList, check for valid data, store in tcl
    forAll(superLoopList, sl)
    {
        timeCluster tc
        (
            readSuperLoopTimes(solverDomain, superLoopList[sl], processor)
        );

        // If there are no time directories, ignore this superLoop
        if (tc.times().size() == 0) continue;

        // Store timeCluster
        tcl.setSize(++nTimeClusters);
        tcl[nTimeClusters - 1] = tc;
    }
    return tcl;
}


Foam::timeClusterList Foam::multiSolver::readAllTimes
(
    const word processor
) const
{
    timeClusterList tcl(0);

    // Loop through solverDomains
    forAll(prefixes_, pf)
    {
        if (prefixes_[pf] == "default") continue;

        timeClusterList tclIn(readSolverDomainTimes(prefixes_[pf], processor));
        tcl.append(tclIn);
    }
    return tcl;
}


bool Foam::multiSolver::loadTimeClusterList
(
    const Foam::timeClusterList& tcl,
    const bool useGlobalTime,
    const bool loadStoreFields
)
{
    if (!nonOverlapping(tcl, useGlobalTime)) return false;

    wordList storeFields;

    forAll(tcl, i)
    {
        fileName currentPath
        (
            findInstancePath(tcl[i], 0).path()
        );

        instantList il(Time::findTimes(currentPath));
        fileName storeFieldsPath
        (
            currentPath/il[Time::findClosestTimeIndex(il, -1.0)].name()
        );

        setSolverDomainPostProcessing(tcl[i].solverDomainName());

        if
        (
            loadStoreFields
         && currentSolverDomainDict_.found("storeFields")
        )
        {
            storeFields = wordList(currentSolverDomainDict_.lookup("storeFields"));
        }
        else
        {
            storeFields.clear();
        }

        forAll(tcl[i].times(), j)
        {
            fileName storeFieldsDestination
            (
                multiDictRegistry_.path()/tcl[i].times()[j].name()
            );

            cp
            (
                currentPath/tcl[i].times()[j].name(),
                multiDictRegistry_.path()
            );
            if (useGlobalTime)
            {
                storeFieldsDestination = multiDictRegistry_.path()/
                    Time::timeName
                (
                    tcl[i].globalValue(j)
                );

                mv
                (
                    multiDictRegistry_.path()/tcl[i].times()[j].name(),
                    storeFieldsDestination
                );
            }
            if
            (
                loadStoreFields
             && (storeFieldsPath != storeFieldsDestination)
            )
            {
                forAll(storeFields, j)
                {
                    cp
                    (
                        storeFieldsPath/storeFields[j],
                        storeFieldsDestination
                    );
                }
            }
        } // end cycle through instants
    } // end cycle through timeClusters
    return true;
}


void Foam::multiSolver::archiveTimeDirs
(
    const Foam::fileName& sourcePath,
    const Foam::fileName& archivePath,
    const Foam::label& purgeWrite
)
{
    if (archivePath.name() == "initial")
    {
        FatalErrorIn("multiSolver::archiveTimeDirs")
            << "Attempting to archive to the 'initial' directory.  This is "
            << "not permitted.  sourcePath = " << sourcePath << ", archivePath"
            << " = " << archivePath
            << abort(FatalError);
    }
    if (exists(archivePath))
    {
        purgeTimeDirs(archivePath);
    }

    mkDir(archivePath);

    // Perform purgeWrite of superLoop directories
    if (purgeWrite)
    {
        labelList allSL(findSuperLoops(archivePath.path()));
        label currentSL(atoi(archivePath.name().c_str()));

        sort(allSL);
        label i = 0;
        while (allSL[i] < currentSL)
        {
            i++;
        }

        for (label j = 1; j <= (i - purgeWrite); j++)
        {
            rmDir(archivePath.path()/name(allSL[j]));
        }
    }

    instantList timeDirs(Time::findTimes(sourcePath));

    forAll(timeDirs, i)
    {
        if (timeDirs[i].name() == "constant") continue;
        mv(sourcePath/timeDirs[i].name(), archivePath/timeDirs[i].name());
    }
}

void Foam::multiSolver::purgeTimeDirs(const Foam::fileName& path)
{
    instantList timeDirs(Time::findTimes(path));

    forAll(timeDirs, i)
    {
        if (timeDirs[i].name() == "constant") continue;
        rmDir(path/timeDirs[i].name());
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
