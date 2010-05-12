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

#include <iostream>
#include "StaticHashTable.H"
#include "IOstreams.H"
#include "IStringStream.H"
#include "OStringStream.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    //for (;;)
    {
    StaticHashTable<double> myTable(10);

    myTable.insert("aaa", 1.0);
    myTable.insert("aba", 2.0);
    myTable.insert("aca", 3.0);
    myTable.insert("ada", 4.0);
    myTable.insert("aeq", 5.0);
    myTable.insert("aaw", 6.0);
    myTable.insert("abs", 7.0);
    myTable.insert("acr", 8.0);
    myTable.insert("adx", 9.0);
    myTable.insert("aec", 10.0);

    Pout<< "Foam output operator:" << nl << endl;
    Pout<< myTable << endl;

    //myTable.erase("aaw");
    //myTable.erase("abs");
    //std::cerr << "Size now:" << myTable.size() << '\n';

    Pout<< "toc:" << nl << endl;
    Pout<< myTable.toc() << endl;


    std::cerr << myTable.find("aaa")() << '\n';
    std::cerr << myTable.find("aba")() << '\n';
    std::cerr << myTable.find("aca")() << '\n';
    std::cerr << myTable.find("ada")() << '\n';
    std::cerr << myTable.find("aeq")() << '\n';
    std::cerr << myTable.find("aaw")() << '\n';
    std::cerr << myTable.find("abs")() << '\n';
    std::cerr << myTable.find("acr")() << '\n';
    std::cerr << myTable.find("adx")() << '\n';
    std::cerr << myTable.find("aec")() << '\n';

    std::cerr << myTable["aaa"] << '\n';

    {
        OStringStream os;

        os << myTable;

        IStringStream is(os.str());

        Pout<< "Foam Istream constructor:" << nl << endl;
        StaticHashTable<double> readTable(is, 100);

        Pout<< readTable << endl;
    }

    std::cerr << "\ncopy construct of table\n" << std::endl;

    StaticHashTable<double> myTable1(myTable);
    Pout<< "myTable1:" << myTable1 << endl;

    std::cerr << "\nassignment of table\n" << std::endl;

    StaticHashTable<double> myTable2(100);
    myTable2.transfer(myTable);

    //Pout<< "myTable:" << myTable << endl;

    forAllConstIter(StaticHashTable<double>, myTable2, iter2)
    {
        std::cerr << *iter2 << '\n';
    }

    std::cerr << "\ntable resize 1\n" << std::endl;

    myTable2.resize(1);
    
    forAllConstIter(StaticHashTable<double>, myTable2, iter2)
    {
        std::cerr << *iter2 << '\n';
    }

    std::cerr << "\ntable size 10000\n" << std::endl;

    myTable2.resize(10000);
    
    forAllConstIter(StaticHashTable<double>, myTable2, iter2)
    {
        std::cerr << *iter2 << '\n';
    }

    }


    std::cerr << "\nBye.\n";

    return 0;
}


// ************************************************************************* //
