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

Description
    Add everything necessary for profiling plus a macro

    Originally proposed in
    http://www.cfd-online.com/Forums/openfoam-bugs/64081-feature-proposal-application-level-profiling.html

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef profiling_H
#define profiling_H

#include "profilingTrigger.H"

// to be used at the beginning of a section to be profiled
// profiling ends automatically at the end of a block
#define addProfile(name) Foam::profilingTrigger profileTriggerFor##name (#name)

// Use this if a description with spaces, colons etc should be added
#define addProfile2(name,descr) Foam::profilingTrigger profileTriggerFor##name (descr)

// this is only needed if profiling should end before the end of a block
#define endProfile(name) profileTriggerFor##name.stop()

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
