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

#include "JobInfo.H"
#include "OSspecific.H"
#include "clock.H"
#include "OFstream.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool JobInfo::writeJobInfo(debug::infoSwitch("writeJobInfo", 0));

JobInfo jobInfo;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
JobInfo::JobInfo()
:
    runningJobsDir_(getEnv("FOAM_JOB_DIR")/"runningJobs"),
    finishedJobsDir_(getEnv("FOAM_JOB_DIR")/"finishedJobs"),
    jobFileName_(hostName() + '.' + Foam::name(pid())),
    runningJobPath_(runningJobsDir_/jobFileName_),
    finishedJobPath_(finishedJobsDir_/jobFileName_)
{
    name() = "JobInfo";

    if (writeJobInfo && Pstream::master())
    {
        if (!dir(runningJobsDir_) && !mkDir(runningJobsDir_))
        {
            FatalErrorIn("JobInfo::JobInfo()")
                << "Cannot make JobInfo directory " << runningJobsDir_
                << Foam::exit(FatalError);
        }

        if (!dir(finishedJobsDir_) && !mkDir(finishedJobsDir_))
        {
            FatalErrorIn("JobInfo::JobInfo()")
                << "Cannot make JobInfo directory " << finishedJobsDir_
                << Foam::exit(FatalError);
        }
    }

    constructed = true;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

JobInfo::~JobInfo()
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        mv(runningJobPath_, finishedJobPath_);
    }

    constructed = false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool JobInfo::write(Ostream& JobInfoFile) const
{
    if (writeJobInfo && Pstream::master())
    {
        if (JobInfoFile.good())
        {
            dictionary::write(JobInfoFile, false);
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return true;
    }
}


void JobInfo::write() const
{
    if (writeJobInfo && Pstream::master())
    {
        if (!write(OFstream(runningJobPath_)()))
        {
            FatalErrorIn("JobInfo::write() const")
                << "Failed to write to JobInfo file "
                << runningJobPath_
                << Foam::exit(FatalError);
        }
    }
}


void JobInfo::end(const word& terminationType)
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        add("cpuTime", cpuTime_.elapsedCpuTime());
        add("endDate", clock::date());
        add("endTime", clock::clockTime());

        if (!found("termination"))
        {
            add("termination", terminationType);
        }

        rm(runningJobPath_);
        write(OFstream(finishedJobPath_)());
    }

    constructed = false;
}


void JobInfo::end()
{
    end("normal");
}


void JobInfo::exit()
{
    end("exit");
}


void JobInfo::abort()
{
    end("abort");
}


void JobInfo::signalEnd() const
{
    if (writeJobInfo && constructed && Pstream::master())
    {
        mv(runningJobPath_, finishedJobPath_);
    }

    constructed = false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
