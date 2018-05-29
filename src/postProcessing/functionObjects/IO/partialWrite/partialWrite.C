/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "partialWrite.H"
#include "dictionary.H"
#include "foamTime.H"
#include "IOobjectList.H"
#include "polyMesh.H"
#include "cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(partialWrite, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::partialWrite::partialWrite
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::partialWrite::~partialWrite()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::partialWrite::read(const dictionary& dict)
{
    dict.lookup("objectNames") >> objectNames_;
    dict.lookup("writeInterval") >> writeInterval_;
    writeInstance_ = 0;

    Info<< type() << " " << name() << ":" << nl
        << "    dumping every " << writeInterval_
        << " th outputTime : " << nl << endl ;
    forAllConstIter(HashSet<word>, objectNames_, iter)
    {
        Info<< ' ' << iter.key();
    }

    if (writeInterval_ < 1)
    {
        FatalIOErrorIn("partialWrite::read(const dictionary&)", dict)
            << "Illegal value for writeInterval " << writeInterval_
            << ". It should be >= 1."
            << exit(FatalIOError);
    }

    // Clear out any previously loaded fields
    vsf_.clear();
    vvf_.clear();
    vSpheretf_.clear();
    vSymmtf_.clear();
    vtf_.clear();

    ssf_.clear();
    svf_.clear();
    sSpheretf_.clear();
    sSymmtf_.clear();
    stf_.clear();

    forAllConstIter(HashSet<word>, objectNames_, iter)
    {
        loadField<scalar>(iter.key(), vsf_, ssf_);
        loadField<vector>(iter.key(), vvf_, svf_);
        loadField<sphericalTensor>(iter.key(), vSpheretf_, sSpheretf_);
        loadField<symmTensor>(iter.key(), vSymmtf_, sSymmtf_);
        loadField<tensor>(iter.key(), vtf_, stf_);
    }
}


void Foam::partialWrite::execute()
{
}


void Foam::partialWrite::end()
{
    //Pout<< "end at time " << obr_.time().timeName() << endl;
    // Do nothing - only valid on write
}


void Foam::partialWrite::timeSet()
{
    if (obr_.time().outputTime())
    {
        writeInstance_++;

        if (writeInstance_ == writeInterval_)
        {
            // Next overall dump corresponds to partial write. Change
            // write options to AUTO_WRITE
            writeInstance_ = 0;

            changeWriteOptions<scalar>(vsf_, ssf_, IOobject::AUTO_WRITE);
            changeWriteOptions<vector>(vvf_, svf_, IOobject::AUTO_WRITE);
            changeWriteOptions<sphericalTensor>
            (
                vSpheretf_,
                sSpheretf_,
                IOobject::AUTO_WRITE
            );
            changeWriteOptions<symmTensor>
            (
                vSymmtf_,
                sSymmtf_,
                IOobject::AUTO_WRITE
            );
            changeWriteOptions<tensor>(vtf_, stf_, IOobject::AUTO_WRITE);
        }
        else
        {
            changeWriteOptions<scalar>(vsf_, ssf_, IOobject::NO_WRITE);
            changeWriteOptions<vector>(vvf_, svf_, IOobject::NO_WRITE);
            changeWriteOptions<sphericalTensor>
            (
                vSpheretf_,
                sSpheretf_,
                IOobject::NO_WRITE
            );
            changeWriteOptions<symmTensor>
            (
                vSymmtf_,
                sSymmtf_,
                IOobject::NO_WRITE
            );
            changeWriteOptions<tensor>(vtf_, stf_, IOobject::NO_WRITE);
        }
    }
}


void Foam::partialWrite::write()
{
    // Do nothing. The fields get written through the
    // standard dump
}


// ************************************************************************* //
