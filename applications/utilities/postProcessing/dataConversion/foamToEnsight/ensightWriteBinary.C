/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "ensightWriteBinary.H"
#include <fstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void writeEnsDataBinary
(
    const char* val,
    std::ofstream& ensFile
)
{
    char buffer[80] = {0};
    strcpy(buffer, val);
    ensFile.write(buffer, 80*sizeof(char));
}


void writeEnsDataBinary
(
    const int val,
    std::ofstream& ensFile
)
{
    ensFile.write(reinterpret_cast<const char*>(&val), sizeof(int));
}


void writeEnsDataBinary
(
    const scalarField& sf,
    std::ofstream& ensightFile
)
{
    if (sf.size())
    {
        List<float> temp(sf.size());

        forAll(sf, i)
        {
            temp[i] = float(sf[i]);
        }

        ensightFile.write
        (
            reinterpret_cast<char*>(temp.begin()),
            sf.size()*sizeof(float)
        );
    }
}


// ************************************************************************* //

