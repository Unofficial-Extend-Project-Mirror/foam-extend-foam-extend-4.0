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
    Macro definitions for declaring ClassName()

\*---------------------------------------------------------------------------*/

#ifndef className_H
#define className_H

#include "word.H"
#include "debug.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Adds typeName information from its argument @a TypeNameString to a class
#define ClassName(TypeNameString)                                             \
    static const char* typeName_() { return TypeNameString; }                 \
    static const ::Foam::word typeName;                                       \
    static int debug

#define NamespaceName(TypeNameString)                                         \
    inline const char* typeName_() { return TypeNameString; }                 \
    extern const ::Foam::word typeName;                                       \
    extern int debug

#define TemplateName(template)                                                \
class template##Name                                                          \
{                                                                             \
public:                                                                       \
                                                                              \
    template##Name()                                                          \
    {}                                                                        \
                                                                              \
    ClassName(#template);                                                     \
};


#define defineTypeNameAndDebug(Type, DebugSwitch)                             \
    const ::Foam::word Type::typeName(Type::typeName_());                     \
    int Type::debug(::Foam::debug::debugSwitch(Type::typeName_(), DebugSwitch));

#ifdef __INTEL_COMPILER
#define defineTemplateTypeNameAndDebugWithName(Type, Name, DebugSwitch)       \
    const ::Foam::word Type::typeName(Name);                                  \
    int Type::debug(::Foam::debug::debugSwitch(Name, DebugSwitch));
#else
#define defineTemplateTypeNameAndDebugWithName(Type, Name, DebugSwitch)       \
    template<>                                                                \
    const ::Foam::word Type::typeName(Name);                                  \
    template<>                                                                \
    int Type::debug(::Foam::debug::debugSwitch(Name, DebugSwitch));
#endif

#define defineTemplateTypeNameAndDebug(Type, DebugSwitch)                     \
    defineTemplateTypeNameAndDebugWithName(Type, #Type, DebugSwitch);

#ifdef __INTEL_COMPILER
#define defineNamedTemplateTypeNameAndDebug(Type, DebugSwitch)                \
    const ::Foam::word Type::typeName(Type::typeName_());                     \
    int Type::debug(::Foam::debug::debugSwitch(Type::typeName_(), DebugSwitch));
#else
#define defineNamedTemplateTypeNameAndDebug(Type, DebugSwitch)                \
    template<>                                                                \
    const ::Foam::word Type::typeName(Type::typeName_());                     \
    template<>                                                                \
    int Type::debug(::Foam::debug::debugSwitch(Type::typeName_(), DebugSwitch));
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
