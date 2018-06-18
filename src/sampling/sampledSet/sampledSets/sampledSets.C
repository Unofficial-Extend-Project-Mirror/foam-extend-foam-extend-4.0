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

#include "sampledSets.H"
#include "dictionary.H"
#include "foamTime.H"
#include "volFields.H"
#include "ListListOps.H"
#include "SortableList.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(sampledSets, 0);
bool sampledSets::verbose_ = false;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSets::combineSampledSets
(
    PtrList<coordSet>& masterSampledSets,
    labelListList& indexSets
)
{
    // Combine sampleSets from processors. Sort by curveDist. Return
    // ordering in indexSets.
    // Note: only master results are valid

    masterSampledSets_.clear();
    masterSampledSets_.setSize(size());
    indexSets_.setSize(size());

    const PtrList<sampledSet>& sampledSets = *this;

    forAll (sampledSets, setI)
    {
        const sampledSet& samplePts = sampledSets[setI];

        // Collect data from all processors
        List<List<point> > gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] = samplePts;
        Pstream::gatherList(gatheredPts);

        List<labelList> gatheredSegments(Pstream::nProcs());
        gatheredSegments[Pstream::myProcNo()] = samplePts.segments();
        Pstream::gatherList(gatheredSegments);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] = samplePts.curveDist();
        Pstream::gatherList(gatheredDist);


        // Combine processor lists into one big list.
        List<point> allPts
        (
            ListListOps::combine<List<point> >
            (
                gatheredPts, accessOp<List<point> >()
            )
        );
        labelList allSegments
        (
            ListListOps::combine<labelList>
            (
                gatheredSegments, accessOp<labelList>()
            )
        );
        scalarList allCurveDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );


        if (Pstream::master() && allCurveDist.size() == 0)
        {
            WarningIn("sampledSets::combineSampledSets(..)")
                << "Sample set " << samplePts.name()
                << " has zero points." << endl;
        }

        // Sort curveDist and use to fill masterSamplePts
        SortableList<scalar> sortedDist(allCurveDist);
        indexSets[setI] = sortedDist.indices();

        masterSampledSets.set
        (
            setI,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                List<point>(UIndirectList<point>(allPts, indexSets[setI])),
                allCurveDist
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::sampledSets
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    PtrList<sampledSet>(),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
    searchEngine_(mesh_),
    interpolationScheme_(word::null),
    writeFormat_(word::null)
{
    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/"postProcessing"/name_;
    }
    else
    {
        outputPath_ = mesh_.time().path()/"postProcessing"/name_;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::~sampledSets()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSets::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


void Foam::sampledSets::execute()
{
    // Do nothing - only valid on write
}


void Foam::sampledSets::end()
{
    // Do nothing - only valid on write
}


void Foam::sampledSets::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::sampledSets::write()
{
    if (size())
    {
        const label nFields = classifyFields();

        if (Pstream::master())
        {
            if (debug)
            {
                Pout<< "timeName = " << mesh_.time().timeName() << nl
                    << "scalarFields    " << scalarFields_ << nl
                    << "vectorFields    " << vectorFields_ << nl
                    << "sphTensorFields " << sphericalTensorFields_ << nl
                    << "symTensorFields " << symmTensorFields_ <<nl
                    << "tensorFields    " << tensorFields_ <<nl;
            }

            if (nFields)
            {
                if (debug)
                {
                    Pout<< "Creating directory "
                        << outputPath_/mesh_.time().timeName()
                            << nl << endl;
                }

                mkDir(outputPath_/mesh_.time().timeName());
            }
        }

        if (nFields)
        {
            sampleAndWrite(scalarFields_);
            sampleAndWrite(vectorFields_);
            sampleAndWrite(sphericalTensorFields_);
            sampleAndWrite(symmTensorFields_);
            sampleAndWrite(tensorFields_);
        }
    }
}


void Foam::sampledSets::read(const dictionary& dict)
{
    dict_ = dict;

    bool setsFound = dict_.found("sets");
    if (setsFound)
    {
        dict_.lookup("fields") >> fieldSelection_;
        clearFieldGroups();

        dict.lookup("interpolationScheme") >> interpolationScheme_;
        dict.lookup("setFormat") >> writeFormat_;

        PtrList<sampledSet> newList
        (
            dict_.lookup("sets"),
            sampledSet::iNew(mesh_, searchEngine_)
        );
        transfer(newList);
        combineSampledSets(masterSampledSets_, indexSets_);

        if (this->size())
        {
            Info<< "Reading set description:" << nl;
            forAll(*this, setI)
            {
                Info<< "    " << operator[](setI).name() << nl;
            }
            Info<< endl;
        }
    }

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample sets:" << nl << "(" << nl;

        forAll(*this, setI)
        {
            Pout<< "  " << operator[](setI) << endl;
        }
        Pout<< ")" << endl;
    }
}


void Foam::sampledSets::correct()
{
    // Reset interpolation
    pointMesh::Delete(mesh_);
    volPointInterpolation::Delete(mesh_);

    searchEngine_.correct();

    PtrList<sampledSet> newList
    (
        dict_.lookup("sets"),
        sampledSet::iNew(mesh_, searchEngine_)
    );
    transfer(newList);
    combineSampledSets(masterSampledSets_, indexSets_);
}


void Foam::sampledSets::updateMesh(const mapPolyMesh&)
{
    correct();
}


void Foam::sampledSets::movePoints(const pointField&)
{
    correct();
}


void Foam::sampledSets::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        correct();
    }
}


// ************************************************************************* //
