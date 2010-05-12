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

#include <iostream>
#include "HashTable.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main()
{
    //for (;;)
    {
    HashTable<double> myTable(100);

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

    myTable.erase("aaw");
    myTable.erase("abs");

    std::cerr << myTable.find("aaa")() << '\n';
    std::cerr << myTable.find("aba")() << '\n';
    std::cerr << myTable.find("aca")() << '\n';
    std::cerr << myTable.find("ada")() << '\n';
    std::cerr << myTable.find("aeq")() << '\n';
    //std::cerr << myTable.find("aaw")() << '\n';
    //std::cerr << myTable.find("abs")() << '\n';
    std::cerr << myTable.find("acr")() << '\n';
    std::cerr << myTable.find("adx")() << '\n';
    std::cerr << myTable.find("aec")() << '\n';

    std::cerr << "\nprint table\n" << std::endl;

    for
    (
        HashTable<double>::iterator iter = myTable.begin();
        iter != myTable.end();
        ++iter
    )
    {
        std::cerr << *iter << '\n';
    }

    std::cerr << "\nprint table\n" << std::endl;

    forAllIter(HashTable<double>, myTable, iter)
    {
        std::cerr << *iter << '\n';
    }

    std::cerr << "\ncopy of table\n" << std::endl;

    HashTable<double> myTable2;
    myTable2 = myTable;

    forAllConstIter(HashTable<double>, myTable2, iter2)
    {
        std::cerr << *iter2 << '\n';
    }

    std::cerr << "\ndelete entries\n" << std::endl;

    forAllIter(HashTable<double>, myTable, iter)
    {
        std::cerr << "deleting " << *iter << '\n';
        myTable.erase(iter);
        std::cerr << "deleted\n";
    }

    forAllConstIter(HashTable<double>, myTable, iter)
    {
        std::cerr << *iter << '\n';
    }
    }

    std::cerr << "\nBye.\n";

    return 0;
}


// ************************************************************************* //
