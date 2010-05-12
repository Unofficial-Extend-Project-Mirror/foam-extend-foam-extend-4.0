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

#include "Time.H"
#include "PstreamReduceOps.H"

#include <sstream>

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(Time, 0);
}


template<>
const char* Foam::NamedEnum<Foam::Time::stopAtControls, 4>::names[] =
{
    "endTime",
    "noWriteNow",
    "writeNow",
    "nextWrite"
};

const Foam::NamedEnum<Foam::Time::stopAtControls, 4>
    Foam::Time::stopAtControlNames_;

template<>
const char* Foam::NamedEnum<Foam::Time::writeControls, 5>::names[] =
{
    "timeStep",
    "runTime",
    "adjustableRunTime",
    "clockTime",
    "cpuTime"
};

const Foam::NamedEnum<Foam::Time::writeControls, 5>
    Foam::Time::writeControlNames_;

Foam::Time::fmtflags Foam::Time::format_(Foam::Time::general);
int Foam::Time::precision_(6);

Foam::word Foam::Time::controlDictName("controlDict");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::Time::adjustDeltaT()
{
    if (writeControl_ == wcAdjustableRunTime)
    {
        scalar timeToNextWrite = max
        (
            0.0,
            (outputTimeIndex_ + 1)*writeInterval_ - (value() - startTime_)
        );

        label nStepsToNextWrite = label(timeToNextWrite/deltaT_ - SMALL) + 1;
        scalar newDeltaT = timeToNextWrite/nStepsToNextWrite;

        // Control the increase of the time step to within a factor of 2
        // and the decrease within a factor of 5.
        if (newDeltaT >= deltaT_)
        {
            deltaT_ = min(newDeltaT, 2.0*deltaT_);
        }
        else
        {
            deltaT_ = max(newDeltaT, 0.2*deltaT_);
        }
    }
}


void Foam::Time::setControls()
{
    // default is to resume calculation from "latestTime"
    word startFrom("latestTime");

    controlDict_.readIfPresent("startFrom", startFrom);

    if (startFrom == "startTime")
    {
        controlDict_.lookup("startTime") >> startTime_;
    }
    else
    {
        // Search directory for valid time directories
        instantList Times = findTimes(path());

        if (startFrom == "firstTime")
        {
            if (Times.size())
            {
                startTime_ = Times[0].value();
            }
        }
        else if (startFrom == "latestTime")
        {
            if (Times.size())
            {
                startTime_ = Times[Times.size()-1].value();
            }
        }
        else
        {
            WarningIn("Time::setControls()")
                << "    expected startTime, firstTime or latestTime"
                << " found '" << startFrom
                << "' in dictionary " << controlDict_.name() << nl
                << "    Setting time to " << startTime_ << endl;
        }
    }

    setTime(startTime_, 0);

    readDict();
    deltaTSave_ = deltaT_;
    deltaT0_ = deltaTSave_;

    if (Pstream::parRun())
    {
        scalar sumStartTime = startTime_;
        reduce(sumStartTime, sumOp<scalar>());
        if
        (
            mag(Pstream::nProcs()*startTime_ - sumStartTime)
          > Pstream::nProcs()*deltaT_/10.0
        )
        {
            FatalErrorIn("Time::setControls()")
                << "Start time is not the same for all processors" << nl
                << "processor " << Pstream::myProcNo() << " has startTime "
                << startTime_ << exit(FatalError);
        }
    }

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (timeDict.readIfPresent("deltaT", deltaTSave_))
    {
        deltaT0_ = deltaTSave_;
    }

    if (timeDict.readIfPresent("index", startTimeIndex_))
    {
        timeIndex_ = startTimeIndex_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Time::Time
(
    const word& controlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    purgeWrite_(0),
    subCycling_(false),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(true),

    readLibs_(controlDict_, "libs"),
    functionObjects_(*this)
{
    setControls();
}


Foam::Time::Time
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        dict
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    purgeWrite_(0),
    subCycling_(false),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(true),

    readLibs_(controlDict_, "libs"),
    functionObjects_(*this)
{
    setControls();
}


Foam::Time::Time
(
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    TimePaths
    (
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    objectRegistry(*this),

    controlDict_
    (
        IOobject
        (
            controlDictName,
            system(),
            *this,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    ),

    startTimeIndex_(0),
    startTime_(0),
    endTime_(0),

    stopAt_(saEndTime),
    writeControl_(wcTimeStep),
    writeInterval_(GREAT),
    purgeWrite_(0),
    subCycling_(false),

    writeFormat_(IOstream::ASCII),
    writeVersion_(IOstream::currentVersion),
    writeCompression_(IOstream::UNCOMPRESSED),
    graphFormat_("raw"),
    runTimeModifiable_(true),

    readLibs_(controlDict_, "libs"),
    functionObjects_(*this)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Time::~Time()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::Time::timeName(const scalar t)
{
    std::ostringstream osBuffer;
    osBuffer.setf(ios_base::fmtflags(format_), ios_base::floatfield);
    osBuffer.precision(precision_);
    osBuffer << t;
    return osBuffer.str();
}


Foam::word Foam::Time::timeName() const
{
    return dimensionedScalar::name();
    //return timeName(timeOutputValue());
}


// Search the construction path for times
Foam::instantList Foam::Time::times() const
{
    return findTimes(path());
}


Foam::word Foam::Time::findInstancePath(const instant& t) const
{
    instantList times = Time::findTimes(path());

    forAllReverse(times, i)
    {
        if (times[i] == t)
        {
            return times[i].name();
        }
    }

    return word::null;
}


Foam::instant Foam::Time::findClosestTime(const scalar t) const
{
    instantList times = Time::findTimes(path());

    // If there is only one time it is "constant" so return it
    if (times.size() == 1)
    {
        return times[0];
    }

    if (t < times[1].value())
    {
        return times[1];
    }
    else if (t > times[times.size() - 1].value())
    {
        return times[times.size() - 1];
    }

    label nearestIndex = -1;
    scalar deltaT = GREAT;

    for (label i=1; i<times.size(); i++)
    {
        scalar diff = mag(times[i].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = i;
        }
    }

    return times[nearestIndex];
}

//
// This should work too,
// if we don't worry about checking "constant" explicitly
//
// Foam::instant Foam::Time::findClosestTime(const scalar t) const
// {
//     instantList times = Time::findTimes(path());
//     label timeIndex = min(findClosestTimeIndex(times, t), 0);
//     return times[timeIndex];
// }

Foam::label Foam::Time::findClosestTimeIndex
(
    const instantList& times,
    const scalar t
)
{
    label nearestIndex = -1;
    scalar deltaT = GREAT;

    forAll (times, i)
    {
        if (times[i].name() == "constant") continue;

        scalar diff = fabs(times[i].value() - t);
        if (diff < deltaT)
        {
            deltaT = diff;
            nearestIndex = i;
        }
    }

    return nearestIndex;
}


Foam::label Foam::Time::startTimeIndex() const
{
    return startTimeIndex_;
}


Foam::dimensionedScalar Foam::Time::startTime() const
{
    return dimensionedScalar("startTime", dimTime, startTime_);
}


Foam::dimensionedScalar Foam::Time::endTime() const
{
    return dimensionedScalar("endTime", dimTime, endTime_);
}


bool Foam::Time::run() const
{
    bool running = value() < (endTime_ - 0.5*deltaT_);

    if (!subCycling_ && !running && timeIndex_ != startTimeIndex_)
    {
        const_cast<functionObjectList&>(functionObjects_).execute();
    }

    return running;
}


bool Foam::Time::end() const
{
    return (value() > (endTime_ + 0.5*deltaT_));
}


void Foam::Time::setTime(const Time& t)
{
    value() = t.value();
    dimensionedScalar::name() = t.dimensionedScalar::name();
    timeIndex_ = t.timeIndex_;
}


void Foam::Time::setTime(const instant& inst, const label newIndex)
{
    value() = inst.value();
    dimensionedScalar::name() = inst.name();
    timeIndex_ = newIndex;

    IOdictionary timeDict
    (
        IOobject
        (
            "time",
            timeName(),
            "uniform",
            *this,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    timeDict.readIfPresent("deltaT", deltaT_);
    timeDict.readIfPresent("deltaT0", deltaT0_);
    timeDict.readIfPresent("index", timeIndex_);
}


void Foam::Time::setTime(const dimensionedScalar& newTime, const label newIndex)
{
    setTime(newTime.value(), newIndex);
}


void Foam::Time::setTime(const scalar newTime, const label newIndex)
{
    value() = newTime;
    dimensionedScalar::name() = timeName(timeToUserTime(newTime));
    timeIndex_ = newIndex;
}


void Foam::Time::setStopAt(const stopAtControls& sa)
{
    stopAt_ = sa;
}


void Foam::Time::setEndTime(const dimensionedScalar& endTime)
{
    setEndTime(endTime.value());
}


void Foam::Time::setEndTime(const scalar endTime)
{
    endTime_ = endTime;
}


void Foam::Time::setDeltaT(const dimensionedScalar& deltaT)
{
    setDeltaT(deltaT.value());
}


void Foam::Time::setDeltaT(const scalar deltaT)
{
    deltaT_ = deltaT;
    deltaTchanged_ = true;
    adjustDeltaT();
}


void Foam::Time::setWriteControl(const writeControls& wc)
{
    writeControl_ = wc;
}


void Foam::Time::setWriteInterval(const scalar writeInterval)
{
    writeInterval_ = writeInterval;
}


Foam::TimeState Foam::Time::subCycle(const label nSubCycles)
{
    subCycling_ = true;

    TimeState ts = *this;
    setTime(*this - deltaT(), (timeIndex() - 1)*nSubCycles);
    deltaT_ /= nSubCycles;
    deltaT0_ /= nSubCycles;
    deltaTSave_ = deltaT0_;

    return ts;
}


void Foam::Time::endSubCycle(const TimeState& ts)
{
    subCycling_ = false;
    TimeState::operator=(ts);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::Time& Foam::Time::operator+=(const dimensionedScalar& deltaT)
{
    return operator+=(deltaT.value());
}


Foam::Time& Foam::Time::operator+=(const scalar deltaT)
{
    readModifiedObjects();

    if (!subCycling_)
    {
        if (timeIndex_ == startTimeIndex_)
        {
            functionObjects_.start();
        }
        else
        {
            functionObjects_.execute();
        }
    }

    setDeltaT(deltaT);
    operator++();

    return *this;
}


Foam::Time& Foam::Time::operator++()
{
    readModifiedObjects();

    if (!subCycling_)
    {
        if (timeIndex_ == startTimeIndex_)
        {
            functionObjects_.start();
        }
        else
        {
            functionObjects_.execute();
        }
    }

    deltaT0_ = deltaTSave_;
    deltaTSave_ = deltaT_;
    setTime(value() + deltaT_, timeIndex_ + 1);

    // If the time is very close to zero reset to zero
    if (mag(value()) < 10*SMALL*deltaT_)
    {
        setTime(0.0, timeIndex_);
    }

    switch(writeControl_)
    {
        case wcTimeStep:
            outputTime_ = !(timeIndex_%label(writeInterval_));
        break;

        case wcRunTime:
        case wcAdjustableRunTime:
        {
            label outputTimeIndex =
                label(((value() - startTime_) + 0.5*deltaT_)/writeInterval_);

            if (outputTimeIndex > outputTimeIndex_)
            {
                outputTime_ = true;
                outputTimeIndex_ = outputTimeIndex;
            }
            else
            {
                outputTime_ = false;
            }
        }
        break;

        case wcCpuTime:
        {
            label outputTimeIndex =
                label(elapsedCpuTime()/writeInterval_);

            if (outputTimeIndex > outputTimeIndex_)
            {
                outputTime_ = true;
                outputTimeIndex_ = outputTimeIndex;
            }
            else
            {
                outputTime_ = false;
            }
        }
        break;

        case wcClockTime:
        {
            label outputTimeIndex = label(elapsedClockTime()/writeInterval_);
            if (outputTimeIndex > outputTimeIndex_)
            {
                outputTime_ = true;
                outputTimeIndex_ = outputTimeIndex;
            }
            else
            {
                outputTime_ = false;
            }
        }
        break;
    };

    if (!end())
    {
        if (stopAt_ == saNoWriteNow)
        {
            endTime_ = value();
        }
        else if (stopAt_ == saWriteNow)
        {
            endTime_ = value();
            outputTime_ = true;
        }
        else if (stopAt_ == saNextWrite && outputTime_ == true)
        {
            endTime_ = value();
        }
    }

    return *this;
}


Foam::Time& Foam::Time::operator++(int)
{
    return operator++();
}


// ************************************************************************* //
