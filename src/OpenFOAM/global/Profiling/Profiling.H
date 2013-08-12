/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright  held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is based on OpenFOAM.

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

Class

Description
    Add everything necessary for Profiling plus a macro

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef Profiling_H
#define Profiling_H

#include "ProfilingTrigger.H"

#define addProfile(name) Foam::ProfilingTrigger profileTriggerFor##name (#name)
// Use this if a description with spaces, colons etc should be added
#define addProfile2(name,descr) Foam::ProfilingTrigger profileTriggerFor##name (descr)
#define endProfile(name) profileTriggerFor##name.stop()

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
