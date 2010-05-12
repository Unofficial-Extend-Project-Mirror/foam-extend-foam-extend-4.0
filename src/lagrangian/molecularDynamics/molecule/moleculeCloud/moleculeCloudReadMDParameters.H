Info<< nl <<  "Reading MD solution parameters:" << endl;

IOdictionary mdSolution
(
    IOobject
    (
        "mdSolution",
        mesh_.time().system(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

integrationMethod_ = integrationMethodNames_.read
(
    mdSolution.lookup("integrationMethod")
);

potentialEnergyLimit_ = readScalar
(
    mdSolution.lookup("potentialEnergyLimit")
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

IOdictionary potentialDict
(
    IOobject
    (
        "potentialDict",
        mesh_.time().system(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

IOdictionary idListDict
(
    IOobject
    (
        "idList",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// ****************************************************************************
// Pair potentials

if (!potentialDict.found("pair"))
{
    FatalErrorIn("moleculeCloudReadPotentials.H") << nl
            << "pair potential specification subDict not found"
            << abort(FatalError);
}

const dictionary& pairDict = potentialDict.subDict("pair");

pairPotentials_.buildPotentials(idListDict, pairDict, mesh_);

if (potentialDict.found("removalOrder"))
{
    List<word> remOrd = potentialDict.lookup("removalOrder");

    removalOrder_.setSize(remOrd.size());

    forAll(removalOrder_, rO)
    {
        removalOrder_[rO] = findIndex(pairPotentials_.idList(), remOrd[rO]);
    }
}

// ****************************************************************************
// Tether potentials

iterator mol(this->begin());

DynamicList<label> tetherIds;

for
(
 mol = this->begin();
 mol != this->end();
 ++mol
)
{
    if (mol().tethered())
    {
        if (findIndex(tetherIds, mol().id()) == -1)
        {
            tetherIds.append
            (
                mol().id()
            );
        }
    }
}

if (Pstream::parRun())
{
    List< labelList > allTetherIds(Pstream::nProcs());

    allTetherIds[Pstream::myProcNo()] = tetherIds;

    Pstream::gatherList(allTetherIds);

    if (Pstream::master())
    {
        DynamicList<label> globalTetherIds;

        forAll(allTetherIds, procN)
        {
            const labelList& procNTetherIds = allTetherIds[procN];

            forAll(procNTetherIds, id)
            {
                if (findIndex(globalTetherIds, procNTetherIds[id]) == -1)
                {
                    globalTetherIds.append
                    (
                        procNTetherIds[id]
                    );
                }
            }
        }

        globalTetherIds.shrink();

        tetherIds = globalTetherIds;
    }

    Pstream::scatter(tetherIds);
}

tetherIds.shrink();

if (tetherIds.size())
{
    if (!potentialDict.found("tether"))
    {
        FatalErrorIn("moleculeCloudReadPotentials.H") << nl
                << "tether potential specification subDict not found"
                << abort(FatalError);
    }

    const dictionary& tetherDict = potentialDict.subDict("tether");

    tetherPotentials_.buildPotentials(idListDict, tetherDict, tetherIds);
}

// ****************************************************************************
// External Forces

gravity_ = vector::zero;

if (potentialDict.found("external"))
{

    Info << nl << "Reading external forces:" << endl;

    const dictionary& externalDict = potentialDict.subDict("external");

    // ************************************************************************
    // gravity

    if (externalDict.found("gravity"))
    {
        gravity_ = externalDict.lookup("gravity");
    }
}

Info << nl << tab << "gravity = " << gravity_ << endl;

