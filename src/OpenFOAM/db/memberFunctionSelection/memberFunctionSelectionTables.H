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

Class
    Foam::memberFunctionSelectionTables

Description
    Macros to enable the easy declaration of member function selection tables.

\*---------------------------------------------------------------------------*/

#include "token.H"

#ifndef memberFunctionSelectionTables_H
#define memberFunctionSelectionTables_H

#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define declareMemberFunctionSelectionTable\
(returnType,baseType,memberFunction,argNames,argList,parList)                 \
                                                                              \
    /* Construct from argList function pointer type */                        \
    typedef returnType (*memberFunction##argNames##MemberFunctionPtr)argList; \
                                                                              \
    /* Construct from argList function table type */                          \
    typedef HashTable                                                         \
        <memberFunction##argNames##MemberFunctionPtr, word, string::hash>     \
        memberFunction##argNames##MemberFunctionTable;                        \
                                                                              \
    /* Construct from argList function pointer table pointer */               \
    static memberFunction##argNames##MemberFunctionTable*                     \
        memberFunction##argNames##MemberFunctionTablePtr_;                    \
                                                                              \
    /* Class to add constructor from argList to table */                      \
    template<class baseType##Type>                                            \
    class add##memberFunction##argNames##MemberFunctionToTable                \
    {                                                                         \
    public:                                                                   \
                                                                              \
        add##memberFunction##argNames##MemberFunctionToTable                  \
        (                                                                     \
            const word& lookup = baseType##Type::typeName                     \
        )                                                                     \
        {                                                                     \
            construct##memberFunction##argNames##MemberFunctionTables();      \
            memberFunction##argNames##MemberFunctionTablePtr_->insert         \
            (                                                                 \
                lookup,                                                       \
                baseType##Type::memberFunction                                \
            );                                                                \
        }                                                                     \
                                                                              \
        ~add##memberFunction##argNames##MemberFunctionToTable()               \
        {                                                                     \
            destroy##memberFunction##argNames##MemberFunctionTables();        \
        }                                                                     \
    };                                                                        \
                                                                              \
    /* Table MemberFunction called from the table add function */             \
    static void construct##memberFunction##argNames##MemberFunctionTables();  \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    static void destroy##memberFunction##argNames##MemberFunctionTables()


#define defineMemberFunctionSelectionTableMemberFunction\
(baseType,memberFunction,argNames)                                            \
                                                                              \
    /* Table MemberFunction called from the table add function */             \
    void baseType::construct##memberFunction##argNames##MemberFunctionTables()\
    {                                                                         \
        static bool constructed = false;                                      \
                                                                              \
        if (!constructed)                                                     \
        {                                                                     \
            baseType::memberFunction##argNames##MemberFunctionTablePtr_       \
                = new baseType::memberFunction##argNames##MemberFunctionTable;\
                                                                              \
            constructed = true;                                               \
        }                                                                     \
    }

#define defineMemberFunctionSelectionTableDestructor\
(baseType,memberFunction,argNames)                                            \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    void baseType::destroy##memberFunction##argNames##MemberFunctionTables()  \
    {                                                                         \
        if (baseType::memberFunction##argNames##MemberFunctionTablePtr_)      \
        {                                                                     \
            delete baseType::memberFunction##argNames##MemberFunctionTablePtr_;\
            baseType::memberFunction##argNames##MemberFunctionTablePtr_ = NULL;\
        }                                                                     \
    }

#define defineMemberFunctionSelectionTablePtr\
(baseType,memberFunction,argNames)                                            \
                                                                              \
    /* Define the constructor function table */                               \
    baseType::memberFunction##argNames##MemberFunctionTable*                  \
        baseType::memberFunction##argNames##MemberFunctionTablePtr_ = NULL

#define defineTemplateMemberFunctionSelectionTablePtr\
(baseType,memberFunction,argNames)                                            \
                                                                              \
    /* Define the constructor function table */                               \
    typename baseType::memberFunction##argNames##MemberFunctionTable*         \
        baseType::memberFunction##argNames##MemberFunctionTablePtr_ = NULL

#define defineMemberFunctionSelectionTable(baseType,memberFunction,argNames)  \
                                                                              \
    defineMemberFunctionSelectionTablePtr                                     \
        (baseType,memberFunction,argNames);                                   \
    defineMemberFunctionSelectionTableMemberFunction                          \
        (baseType,memberFunction,argNames)                                    \
    defineMemberFunctionSelectionTableDestructor                              \
        (baseType,memberFunction,argNames)

#define defineTemplateMemberFunctionSelectionTable\
(baseType,memberFunction,argNames)                                            \
                                                                              \
    template<>                                                                \
    defineMemberFunctionSelectionTablePtr                                     \
        (baseType,memberFunction,argNames);                                   \
    template<>                                                                \
    defineMemberFunctionSelectionTableMemberFunction                          \
        (baseType,memberFunction,argNames)                                    \
    template<>                                                                \
    defineMemberFunctionSelectionTableDestructor                              \
        (baseType,memberFunction,argNames)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
