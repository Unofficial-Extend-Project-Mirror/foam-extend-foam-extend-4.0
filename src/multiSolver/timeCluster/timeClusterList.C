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

#include "timeClusterList.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::timeClusterList::typeName = "timeClusterList";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeClusterList::timeClusterList()
:
    List<timeCluster>()
{}


Foam::timeClusterList::timeClusterList(const Foam::label size)
:
    List<timeCluster>(size)
{}


Foam::timeClusterList::timeClusterList(const timeCluster& tcIn)
:
    List<timeCluster>(1, tcIn)
{}


Foam::timeClusterList::timeClusterList(const Foam::label size, const Foam::timeCluster& tcIn)
:
    List<timeCluster>(size, tcIn)
{}


Foam::timeClusterList::timeClusterList(const labelList& subIndices, const timeClusterList& tclIn)
:
    List<timeCluster>(subIndices.size())
{
    forAll(subIndices, i)
    {
        this->operator[](i) = tclIn[subIndices[i]];
    }
}


Foam::timeClusterList::timeClusterList(Foam::Istream& is)
:
    List<timeCluster>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeClusterList::globalSort()
{
    Foam::sort(*this, timeCluster::less());
}

void Foam::timeClusterList::append(const timeClusterList& tclIn)
{
    label wasSize = this->size();
    this->setSize(tclIn.size() + this->size());
    
    for (label i = 0; i < tclIn.size(); i++)
    {
        this->operator[](i + wasSize) = tclIn[i];
    }
}


void Foam::timeClusterList::append(const timeCluster& tcIn)
{
    label wasSize = this->size();
    this->setSize(this->size() + 1);
    
    this->operator[](wasSize) = tcIn;
}


bool Foam::timeClusterList::purgeEmpties()
{
    if (!this->size()) return false;
    
    label empties(0);
    for (label i = 0; i < this->size(); i++)
    {
        if (!this->operator[](i).times().size())
        {
            empties++;
            continue;
        }
        if (empties)
        {
            this->operator[](i - empties) = this->operator[](i);
        }
    }
    if (empties)
    {
        this->setSize(this->size() - empties);
    }
    if (!this->size()) return false;
    return true;
}

Foam::timeClusterList Foam::timeClusterList::selectiveSubList
(
    const labelList& indices
) const
{
    timeClusterList tcl(indices.size());
    
    forAll(indices, i)
    {
        if (indices[i] > this->size())
        {
            FatalErrorIn("timeClusterList::selectiveSubList")
                << "Out of range index passed to this function.  Indices "
                << "passed are: \n" << indices << "\nFailure at index " << i
                << ", with value " << indices[i] << ".\n This timeClusterList "
                << "has size " << this->size() << "."
                << abort(FatalError);
        }

        tcl[i] = this->operator[](indices[i]);
    }
    return tcl;
}



// ************************************************************************* //
