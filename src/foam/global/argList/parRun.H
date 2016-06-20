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

Class
    Foam::ParRunControl

Description
    Helper class for initializing parallel jobs from the command arguments.

\*---------------------------------------------------------------------------*/

#ifndef parRun_H
#define parRun_H

#include "OPstream.H"
#include "IPstream.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ParRunControl Declaration
\*---------------------------------------------------------------------------*/

class ParRunControl
{
    bool RunPar;

public:

    ParRunControl()
    :
        RunPar(false)
    {}

    ~ParRunControl()
    {
        if (RunPar)
        {
            Info<< "Finalising parallel run" << endl;
            Pstream::exit(0);
        }
    }

    void runPar(int& argc, char**& argv)
    {
        RunPar = true;

        if (!Pstream::init(argc, argv))
        {
            Info<< "Failed to start parallel run" << endl;
            Pstream::exit(1);
        }
    }

    bool parRun() const
    {
        return RunPar;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
