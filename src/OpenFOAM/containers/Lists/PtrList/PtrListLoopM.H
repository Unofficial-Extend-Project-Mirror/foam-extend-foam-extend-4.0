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
    simple generic PtrList MACROS for looping

\*---------------------------------------------------------------------------*/

#ifndef PtrListLoop_H
#define PtrListLoop_H

#include "undefListLoopM.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Element access looping using [] for vector and parallel machines

#define List_FOR_ALL(f, i)      \
        forAll(f, i)            \
        {                       \

#define List_END_FOR_ALL  }

#define List_ELEM(f, fp, i) ((f)[i])

#define List_ACCESS(type, f, fp)
#define List_CONST_ACCESS(type, f, fp)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
