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

Description

\*---------------------------------------------------------------------------*/

#include "cellTable.H"
#include "IOMap.H"
#include "OFstream.H"
#include "wordList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::cellTable::defaultMaterial_ = "fluid";

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::cellTable::zoneMap() const
{
    label maxId = 0;
    forAllConstIter(Map<dictionary>, *this, iter)
    {
        if (maxId < iter.key())
        {
            maxId = iter.key();
        }
    }

    label zoneI = 0;
    labelList list(maxId+1, -1);

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        list[iter.key()] = zoneI++;
    }

    return list;
}


Foam::wordList Foam::cellTable::namesList() const
{
    Map<word> lookup = names();
    wordList lst(lookup.size());

    label n = 0;
    forAllConstIter(Map<word>, lookup, iter)
    {
        lst[n] = iter();
    }

    lst.setSize(n);

    return lst;
}


void Foam::cellTable::addDefaults()
{
    forAllIter(Map<dictionary>, *this, iter)
    {
        if (!iter().found("MaterialType"))
        {
            iter().add("MaterialType", defaultMaterial_);
        }
    }
}


void Foam::cellTable::setEntry
(
    const label& id,
    const word& keyWord,
    const word& value
)
{
    dictionary dict;
    dict.add(keyWord, value);

    iterator iter = find(id);
    if (iter != end())
    {
        iter().merge(dict);
    }
    else
    {
        insert(id, dict);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellTable::cellTable()
:
    Map<dictionary>()
{}


Foam::cellTable::cellTable
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
:
    Map<dictionary>()
{
    readDict(registry, name, instance);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellTable::~cellTable()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::cellTable::append(const dictionary& dict)
{
    label maxId = -1;
    forAllConstIter(Map<dictionary>, *this, iter)
    {
        if (maxId < iter.key())
        {
            maxId = iter.key();
        }
    }

    insert(++maxId, dict);
    return maxId;
}


Foam::Map<Foam::word> Foam::cellTable::names() const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word theName = "cellTable_" + Foam::name(iter.key());

        iter().readIfPresent("Label", theName);

        lookup.insert(iter.key(), theName);
    }

    return lookup;
}


Foam::word Foam::cellTable::name(const label& id) const
{
    const_iterator iter = find(id);

    word theName = "cellTable_" + Foam::name(id);
    if (iter != end())
    {
        iter().readIfPresent("Label", theName);
    }

    return theName;
}


Foam::label Foam::cellTable::findIndex(const word& name) const
{
    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word theName;
        if (iter().readIfPresent("Label", theName))
        {
            if (theName == name)
            {
                return iter.key();
            }
        }
    }

    return -1;
}


Foam::Map<Foam::word> Foam::cellTable::selectType
(
    const word& materialType
) const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word matl(defaultMaterial_);

        iter().readIfPresent("MaterialType", matl);

        if (matl == materialType)
        {
            word theName = "cellTable_" + Foam::name(iter.key());

            iter().readIfPresent("Label", theName);

            lookup.insert(iter.key(), theName);
        }
    }

    return lookup;
}


Foam::Map<Foam::word> Foam::cellTable::fluids() const
{
    return selectType("fluid");
}


Foam::Map<Foam::word> Foam::cellTable::solids() const
{
    return selectType("solid");
}


Foam::Map<Foam::word> Foam::cellTable::shells() const
{
    return selectType("shell");
}


Foam::Map<Foam::word> Foam::cellTable::materialTypes() const
{
    Map<word> lookup;

    forAllConstIter(Map<dictionary>, *this, iter)
    {
        word matlType(defaultMaterial_);

        iter().readIfPresent("MaterialType", matlType);

        lookup.insert(iter.key(), matlType);
    }

    return lookup;
}


void Foam::cellTable::setMaterial(const label& id, const word& matlType)
{
    setEntry(id, "MaterialType", matlType);
}


void Foam::cellTable::setName(const label& id, const word& name)
{
    setEntry(id, "Label", name);
}


void Foam::cellTable::setName(const label& id)
{
    iterator iter = find(id);

    if (iter == end() || !iter().found("Label"))
    {
        setName(id, "cellTable_" + ::Foam::name(id));
    }
}


void Foam::cellTable::readDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
{
    clear();

    // read constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (ioObj.headerOk())
    {
        *this = ioObj;
        addDefaults();
    }
    else
    {
        Info<< "no constant/cellTable information available" << endl;
    }
}


void Foam::cellTable::writeDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
) const
{
    // write constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    ioObj.note() = "persistent data for thirdParty mesh <-> OpenFOAM translation";

    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);
    os << *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cellTable::operator=(const cellTable& rhs)
{
    Map<dictionary>::operator=(rhs);
    addDefaults();
}


void Foam::cellTable::operator=(const Map<dictionary>& rhs)
{
    Map<dictionary>::operator=(rhs);
    addDefaults();
}


void Foam::cellTable::operator=(const polyMesh& mesh)
{
    Map<dictionary> zoneDict;

    // create cellTableId and cellTable based on cellZones
    label nZoneCells = 0;

    wordList zoneNames = mesh.cellZones().names();
    label unZonedType = zoneNames.size() + 1;

    // do cell zones
    forAll(mesh.cellZones(), zoneI)
    {
        const cellZone& cZone = mesh.cellZones()[zoneI];
        nZoneCells += cZone.size();

        dictionary dict;
        dict.add("Label", zoneNames[zoneI]);
        zoneDict.insert(zoneI + 1, dict);
    }

    // collect unzoned cells
    // special case: no zones at all - do entire mesh
    if (nZoneCells == 0)
    {
        zoneDict.clear();
        unZonedType = 1;
    }

    if (mesh.nCells() > nZoneCells)
    {
        zoneDict.insert
        (
            unZonedType,
            dictionary
            (
                IStringStream("Label cells;")()
            )
        );
    }

    Map<dictionary>::operator=(zoneDict);
    addDefaults();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

void Foam::cellTable::addCellZones
(
    polyMesh& mesh,
    const labelList& tableIds
) const
{
    labelList typeToZone = zoneMap();
    wordList  zoneNames = namesList();

    List<DynamicList<label> > zoneCells(size());

    forAll(tableIds, cellI)
    {
        label zoneI = typeToZone[tableIds[cellI]];
        if (zoneI >= 0)
        {
            zoneCells[zoneI].append(cellI);
        }
    }

    // avoid empty zones
    labelList zoneUsed(zoneCells.size());

    label nZone = 0;
    forAll(zoneCells, zoneI)
    {
        zoneCells[zoneI].shrink();
        if (zoneCells[zoneI].size() > 0)
        {
            zoneUsed[nZone++] = zoneI;
        }
    }
    zoneUsed.setSize(nZone);

    mesh.cellZones().clear();
    if (nZone <= 1)
    {
        Info<< "cellZones not used" << endl;
        return;
    }
    mesh.cellZones().setSize(nZone);

    forAll(zoneUsed, zoneI)
    {
        const label origZoneI = zoneUsed[zoneI];

        Info<< "cellZone " << zoneI
            << " (size: " << zoneCells[origZoneI].size() << ") name: "
            << zoneNames[origZoneI] << endl;

        mesh.cellZones().set
        (
            zoneI,
            new cellZone
            (
                zoneNames[origZoneI],
                zoneCells[origZoneI],
                zoneI,
                mesh.cellZones()
            )
        );
    }
    mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;
}


void Foam::cellTable::combine(const dictionary& dict, labelList& tableIds)
{
    if (!dict.size())
    {
        return;
    }

    bool remap = false;
    labelList mapping(identity(max(this->toc()) + 1));

    forAllConstIter (dictionary, dict, iter)
    {
        wordList  zoneNames(iter().stream());
        labelList zoneIndex(zoneNames.size());

        label nElem = 0;
        forAll (zoneNames, zoneI)
        {
            zoneIndex[nElem] = this->findIndex(zoneNames[zoneI]);
            if (zoneIndex[nElem] >= 0)
            {
                if (zoneI != nElem)
                {
                    zoneNames[nElem] = zoneNames[zoneI];
                }
                ++nElem;
            }
        }

        zoneIndex.setSize(nElem);
        zoneNames.setSize(nElem);

        if (nElem)
        {
            remap = true;
            label targetId = this->findIndex(iter().keyword());

            Info<< "combine cellTable: " << iter().keyword();
            if (targetId >= 0)
            {
                Info<< " += (";
            }
            else
            {
                Info<< " = (";
            }
            forAll (zoneNames, zoneI)
            {
                Info<< " " << zoneNames[zoneI];
            }
            Info<< " )" << endl;

            // re-use the first element if possible
            if (targetId < 0)
            {
                targetId = min(zoneIndex);
                dictionary newDict(operator[](targetId));

                newDict.remove("Label");
                newDict.add("Label", iter().keyword());
                this->set(targetId, newDict);
            }

            forAll (zoneIndex, zoneI)
            {
                label idx = zoneIndex[zoneI];
                if (idx != targetId && idx >= 0)
                {
                    mapping[idx] = targetId;
                    this->erase(idx);
                }
            }
        }
    }

    if (remap)
    {
        inplaceRenumber(mapping, tableIds);
    }
}


// ************************************************************************* //
