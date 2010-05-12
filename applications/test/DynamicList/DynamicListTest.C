/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "DynamicList.H"
#include "IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<DynamicList<label, 1, 0> > ldl(2);

    ldl[0](0) = 0;
    ldl[0](2) = 2;
    ldl[0](3) = 3;
    ldl[0](1) = 1;

    ldl[0].setSize(5);     // increase allocated size
    ldl[1].setSize(10);    // increase allocated size
    ldl[1](2) = 2;

    ldl[1] = 3;

    Info<< "<ldl>" << ldl << "</ldl>" << nl << "sizes: ";
    forAll(ldl, i)
    {
        Info<< " " << ldl[i].size() << "/" << ldl[i].allocSize();
    }
    Info<< endl;

    List<List<label> > ll(2);
    ll[0].transfer(ldl[0]);
    ll[1].transfer(ldl[1].shrink());

    Info<< "<ldl>" << ldl << "</ldl>" << nl << "sizes: ";
    forAll(ldl, i)
    {
        Info<< " " << ldl[i].size() << "/" << ldl[i].allocSize();
    }
    Info<< endl;

    Info<< "<ll>" << ll << "</ll>" << nl << endl;


    // test the transfer between DynamicLists
    DynamicList<label, 1, 0> dlA;
    DynamicList<label, 1, 0> dlB;

    for (label i = 0; i < 5; i++)
    {
        dlA.append(i);
    }
    dlA.setSize(10);

    Info<< "<dlA>" << dlA << "</dlA>" << nl << "sizes: "
        << " " << dlA.size() << "/" << dlA.allocSize() << endl;

    dlB.transfer(dlA);

    // provokes memory error if previous transfer did not maintain
    // the correct allocated space

    // currently fails in OpenFOAM-1.5.x
    /*
     * dlB[6] = 6;
     */

    Info<< "Transferred to dlB" << endl;
    Info<< "<dlA>" << dlA << "</dlA>" << nl << "sizes: "
        << " " << dlA.size() << "/" << dlA.allocSize() << endl;
    Info<< "<dlB>" << dlB << "</dlB>" << nl << "sizes: "
        << " " << dlB.size() << "/" << dlB.allocSize() << endl;


    return 0;
}


// ************************************************************************* //
