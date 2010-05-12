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

InClass
    Foam::memberFunctionSelectionTables

Description
    Macros to enable the easy insertion into member function selection tables.

\*---------------------------------------------------------------------------*/

#ifndef addToMemberFunctionSelectionTable_H
#define addToMemberFunctionSelectionTable_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addToMemberFunctionSelectionTable\
(baseType,thisType,memberFunction,argNames)                                   \
                                                                              \
    /* Add the thisType constructor function to the table */                  \
    baseType::add##memberFunction##argNames##MemberFunctionToTable<thisType>  \
    add##thisType##memberFunction##argNames##MemberFunctionTo##baseType##Table_

#define addNamedToMemberFunctionSelectionTable\
(baseType,thisType,memberFunction,argNames,lookup)                            \
                                                                              \
    /* Add the thisType constructor function to the table */                  \
    baseType::add##memberFunction##argNames##MemberFunctionToTable<thisType>  \
    add##thisType##memberFunction##argNames##MemberFunctionTo##baseType##Table_(#lookup)

#define addTemplateToMemberFunctionSelectionTable\
(baseType,thisType,Targ,memberFunction,argNames)                              \
                                                                              \
    /* Add the thisType constructor function to the table */                  \
    baseType::add##memberFunction##argNames##MemberFunctionToTable            \
    <thisType<Targ> >                                                         \
    add##thisType##Targ##memberFunction##argNames##MemberFunctionTo##baseType##Table_


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
