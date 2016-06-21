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

#include "LUscalarMatrix.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::LUscalarMatrix::solve(Field<Type>& sourceSol) const
{
    if (Pstream::parRun())
    {
        Field<Type> completeSourceSol(n());

        if (Pstream::master())
        {
            typename Field<Type>::subField
            (
                completeSourceSol,
                sourceSol.size()
            ).assign(sourceSol);

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream::read
                (
                    Pstream::scheduled,
                    slave,
                    reinterpret_cast<char*>
                    (
                        &(completeSourceSol[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave + 1] - procOffsets_[slave])*sizeof(Type)
                );
            }
        }
        else
        {
            OPstream::write
            (
                Pstream::scheduled,
                Pstream::masterNo(),
                reinterpret_cast<const char*>(sourceSol.begin()),
                sourceSol.byteSize()
            );
        }

        if (Pstream::master())
        {
            LUBacksubstitute(*this, pivotIndices_, completeSourceSol);

            sourceSol = typename Field<Type>::subField
            (
                completeSourceSol,
                sourceSol.size()
            );

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                OPstream::write
                (
                    Pstream::blocking,
                    slave,
                    reinterpret_cast<const char*>
                    (
                        &(completeSourceSol[procOffsets_[slave]])
                    ),
                    (procOffsets_[slave + 1] - procOffsets_[slave])*sizeof(Type)
                );
            }
        }
        else
        {
            IPstream::read
            (
                Pstream::blocking,
                Pstream::masterNo(),
                reinterpret_cast<char*>(sourceSol.begin()),
                sourceSol.byteSize()
            );
        }
    }
    else
    {
        LUBacksubstitute(*this, pivotIndices_, sourceSol);
    }
}


// ************************************************************************* //
