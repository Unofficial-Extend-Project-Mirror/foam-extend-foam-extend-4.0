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

Description
    Macro definitions for declaring ClassName(), NamespaceName(), etc.

\*---------------------------------------------------------------------------*/

#ifndef className_H
#define className_H

#include "word.H"
#include "debug.H"
#include "debugSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declarations (without debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Add typeName information from argument @a TypeNameString to a class.
//  Without debug information
#define ClassNameNoDebug(TypeNameString)                                      \
    static const char* typeName_() { return TypeNameString; }                 \
    static const ::Foam::word typeName

//- Add typeName information from argument @a TypeNameString to a namespace.
//  Without debug information.
#define NamespaceNameNoDebug(TypeNameString)                                  \
    inline const char* typeName_() { return TypeNameString; }                 \
    extern const ::Foam::word typeName

//- Add typeName information from argument @a TemplateNameString to a template class.
//  Without debug information.
#define TemplateNameNoDebug(TemplateNameString)                               \
class TemplateNameString##Name                                                \
{                                                                             \
public:                                                                       \
    TemplateNameString##Name() {}                                             \
    ClassNameNoDebug(#TemplateNameString);                                    \
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// declarations (with debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Add typeName information from argument @a TypeNameString to a class.
//  Also declares debug information.
#define ClassName(TypeNameString)                                             \
    ClassNameNoDebug(TypeNameString);                                         \
    static Foam::debug::debugSwitch debug;

//- Add typeName information from argument @a TypeNameString to a namespace.
//  Also declares debug information.
#define NamespaceName(TypeNameString)                                         \
    NamespaceNameNoDebug(TypeNameString);                                     \
    extern Foam::debug::debugSwitch debug;

//- Add typeName information from argument @a TypeNameString to a template class.
//  Also declares debug information.
#define TemplateName(TemplateNameString)                                      \
class TemplateNameString##Name                                                \
{                                                                             \
public:                                                                       \
    TemplateNameString##Name() {}                                             \
    ClassName(#TemplateNameString);                                           \
}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// definitions (without debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Define the typeName, with alternative lookup as @a Name
#define defineTypeNameWithName(Type, Name)                                    \
    const ::Foam::word Type::typeName(Name)

//- Define the typeName
#define defineTypeName(Type)                                                  \
    defineTypeNameWithName(Type, Type::typeName_())

#ifdef __INTEL_COMPILER
//- Define the typeName as @a Name for template classes
# define defineTemplateTypeNameWithName(Type, Name)                           \
    defineTypeNameWithName(Type, Name)
#else
//- Define the typeName as @a Name for template classes
# define defineTemplateTypeNameWithName(Type, Name)                           \
    template<>                                                                \
    defineTypeNameWithName(Type, Name)
#endif

//- Define the typeName for template classes, useful with typedefs
#define defineTemplateTypeName(Type)                                          \
    defineTemplateTypeNameWithName(Type, #Type)

//- Define the typeName directly for template classes
#define defineNamedTemplateTypeName(Type)                                     \
    defineTemplateTypeNameWithName(Type, Type::typeName_())



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// definitions (debug information only)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Define the debug information, lookup as @a Name
#define defineDebugSwitchWithName(Type, Name, DebugSwitch, SwitchDescr)       \
    ::Foam::debug::debugSwitch                                                \
    Type::debug(std::string(Name), DebugSwitch, SwitchDescr)

//- Define the debug information
#define defineDebugSwitch(Type, DebugSwitch, SwitchDescr)                     \
    defineDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch, SwitchDescr);

#ifdef __INTEL_COMPILER
//- Define the debug information for templates, lookup as @a Name
# define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)           \
    defineDebugSwitchWithName(Type, Name, DebugSwitch, "");

#else
//- Define the debug information for templates, lookup as @a Name
# define defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)           \
    template<>                                                                \
    defineDebugSwitchWithName(Type, Name, DebugSwitch, "");
#endif

//- Define the debug information for templates
//  Useful with typedefs
#define defineTemplateDebugSwitch(Type, DebugSwitch)                          \
    defineTemplateDebugSwitchWithName(Type, #Type, DebugSwitch)

//- Define the debug information directly for templates
#define defineNamedTemplateDebugSwitch(Type, DebugSwitch)                     \
    defineTemplateDebugSwitchWithName(Type, Type::typeName_(), DebugSwitch)



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// definitions (with debug information)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//- Define the typeName and debug information
#define defineTypeNameAndDebug(Type, DebugSwitch)                             \
    defineTypeName(Type);                                                     \
    defineDebugSwitch(Type, DebugSwitch, "")

//- Define the typeName and debug information + description
#define defineTypeNameAndDebugWithDescription(Type, DebugSwitch, SwitchDescr) \
    defineTypeName(Type);                                                     \
    defineDebugSwitch(Type, DebugSwitch, SwitchDescr)

//- Define the typeName and debug information, lookup as @a Name
#define defineTemplateTypeNameAndDebugWithName(Type, Name, DebugSwitch)       \
    defineTemplateTypeNameWithName(Type, Name);                               \
    defineTemplateDebugSwitchWithName(Type, Name, DebugSwitch)

//- Define the typeName and debug information for templates, useful with typedefs
#define defineTemplateTypeNameAndDebug(Type, DebugSwitch)                     \
    defineTemplateTypeNameAndDebugWithName(Type, #Type, DebugSwitch)

//- Define the typeName and debug information for templates
#define defineNamedTemplateTypeNameAndDebug(Type, DebugSwitch)                \
    defineNamedTemplateTypeName(Type);                                        \
    defineNamedTemplateDebugSwitch(Type, DebugSwitch)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
