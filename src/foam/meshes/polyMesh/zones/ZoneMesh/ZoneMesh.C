/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "ZoneMesh.H"
#include "entry.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType, class MeshType>
void ZoneMesh<ZoneType, MeshType>::calcZoneMap() const
{
    // It is an error to attempt to recalculate zoneMap
    // if the pointer is already set
    if (zoneMapPtr_)
    {
        FatalErrorIn("void ZoneMesh<ZoneType>::calcZoneMap() const")
            << "zone map already calculated"
            << abort(FatalError);
    }
    else
    {
        // Count number of objects in all zones
        label nObjects = 0;

        forAll (*this, zoneI)
        {
            nObjects += this->operator[](zoneI).size();
        }

        zoneMapPtr_ = new Map<label>(2*nObjects);
        Map<label>& zm = *zoneMapPtr_;

        // Fill in objects of all zones into the map.  The key is the global
        // object index and the result is the zone index
        forAll (*this, zoneI)
        {
            const labelList& zoneObjects = this->operator[](zoneI);

            forAll (zoneObjects, objI)
            {
                zm.insert(zoneObjects[objI], zoneI);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Read constructor given IOobject and a MeshType reference
template<class ZoneType, class MeshType>
ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(NULL)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        PtrList<ZoneType>& zones = *this;

        // Read zones
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        zones.setSize(patchEntries.size());

        forAll(zones, zoneI)
        {
            zones.set
            (
                zoneI,
                ZoneType::New
                (
                    patchEntries[zoneI].keyword(),
                    patchEntries[zoneI].dict(),
                    zoneI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "ZoneMesh::ZoneMesh"
            "(const IOobject&, const MeshType&)"
        );

        close();
    }
    else
    {
        // No files found.  Force a write of zero-sized zones
        // write();
    }
}


// Construct given size. Zones will be set later
template<class ZoneType, class MeshType>
ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh,
    const label size
)
:
    PtrList<ZoneType>(size),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
ZoneMesh<ZoneType, MeshType>::~ZoneMesh()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map of zones for quick zone lookup
template<class ZoneType, class MeshType>
const Map<label>& ZoneMesh<ZoneType, MeshType>::zoneMap() const
{
    if (!zoneMapPtr_)
    {
        calcZoneMap();
    }

    return *zoneMapPtr_;
}


// Given a global object index, return the zone it is in.
// If object does not belong to any zones, return -1
template<class ZoneType, class MeshType>
label ZoneMesh<ZoneType, MeshType>::whichZone(const label objectIndex) const
{
    const Map<label>& zm = zoneMap();
    Map<label>::const_iterator zmIter = zm.find(objectIndex);

    if (zmIter == zm.end())
    {
        return -1;
    }
    else
    {
        return zmIter();
    }
}


// Return a list of zone names
template<class ZoneType, class MeshType>
wordList ZoneMesh<ZoneType, MeshType>::types() const
{
    const PtrList<ZoneType>& zones = *this;

    wordList t(zones.size());

    forAll (zones, zoneI)
    {
        t[zoneI] = zones[zoneI].type();
    }

    return t;
}


// Return a list of zone names
template<class ZoneType, class MeshType>
wordList ZoneMesh<ZoneType, MeshType>::names() const
{
    const PtrList<ZoneType>& zones = *this;

    wordList t(zones.size());

    forAll (zones, zoneI)
    {
        t[zoneI] = zones[zoneI].name();
    }

    return t;
}


template<class ZoneType, class MeshType>
label ZoneMesh<ZoneType, MeshType>::findZoneID(const word& zoneName) const
{
    const PtrList<ZoneType>& zones = *this;

    forAll (zones, zoneI)
    {
        if (zones[zoneI].name() == zoneName)
        {
            return zoneI;
        }
    }

    // Zone not found
    if (debug)
    {
        InfoIn
        (
            "label ZoneMesh<ZoneType>::findZoneID(const word& zoneName) const"
        )   << "Zone named " << zoneName << " not found.  "
            << "List of available zone names: " << names() << endl;
    }

    // Return for zone not found
    return -1;
}


template<class ZoneType, class MeshType>
void ZoneMesh<ZoneType, MeshType>::clearAddressing()
{
    deleteDemandDrivenData(zoneMapPtr_);

    PtrList<ZoneType>& zones = *this;

    forAll (zones, zoneI)
    {
        zones[zoneI].clearAddressing();
    }
}


template<class ZoneType, class MeshType>
void ZoneMesh<ZoneType, MeshType>::clear()
{
    clearAddressing();
    PtrList<ZoneType>::clear();
}


// Check zone definition
template<class ZoneType, class MeshType>
bool ZoneMesh<ZoneType, MeshType>::checkDefinition(const bool report) const
{
    bool inError = false;

    const PtrList<ZoneType>& zones = *this;

    forAll (zones, zoneI)
    {
        inError |= zones[zoneI].checkDefinition(report);
    }
    return inError;
}


// Correct zone mesh after moving points
template<class ZoneType, class MeshType>
void ZoneMesh<ZoneType, MeshType>::movePoints(const pointField& p)
{
    PtrList<ZoneType>& zones = *this;

    forAll (zones, zoneI)
    {
        zones[zoneI].movePoints(p);
    }
}


template<class ZoneType, class MeshType>
void ZoneMesh<ZoneType, MeshType>::updateMesh()
{
    this->clearAddressing();

    PtrList<ZoneType>& zones = *this;

    forAll (zones, zoneI)
    {
        zones[zoneI].updateMesh();
    }
}


// writeData member function required by regIOobject
template<class ZoneType, class MeshType>
bool ZoneMesh<ZoneType, MeshType>::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Ostream& operator<<(Ostream& os, const ZoneMesh<ZoneType, MeshType>& zones)
{
    os  << zones.size() << nl << token::BEGIN_LIST;

    forAll(zones, zoneI)
    {
        zones[zoneI].writeDict(os);
    }

    os  << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
