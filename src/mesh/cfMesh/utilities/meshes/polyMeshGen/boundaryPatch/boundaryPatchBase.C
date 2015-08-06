/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "boundaryPatchBase.H"
#include "Ostream.H"
#include "Istream.H"
#include "token.H"
#include "IOPtrList.H"
#include "dictionary.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebugWithName
(
    IOPtrList<boundaryPatchBase>,
    "polyBoundaryMesh",
    0
);
    
defineTypeNameAndDebug(boundaryPatchBase, 0);
defineRunTimeSelectionTable(boundaryPatchBase, dictionary);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<boundaryPatchBase> boundaryPatchBase::New
(
    const word& name,
    const dictionary& dict
)
{
    word type(dict.lookup("type"));
    //- check the type of processor. Allowed types are processor and patch
    //- Other patch types are treated as ordinary patches
    if( type != "processor" )
        type = "patch";
    
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if( cstrIter == dictionaryConstructorTablePtr_->end() )
    {
        FatalIOErrorIn
        (
            "boundaryPatchBase::New(const word&, const dictionary&)",
            dict
        )   << "Unknown boundaryPatchBase type " << type << nl << nl
            << "Valid boundaryPatchBase types are :" << nl
            << "[default: " << typeName_() << "]"
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }
    
    return autoPtr<boundaryPatchBase>(cstrIter()(name, dict));
}

autoPtr<boundaryPatchBase> boundaryPatchBase::New
(
    Istream& is
)
{
    word name(is);
    dictionary dict(is);

    return boundaryPatchBase::New(name, dict);
}
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

boundaryPatchBase::boundaryPatchBase
(
    const word& n,
    const word& t,
    const label nF,
    const label sF
)
:
    name_(n),
    type_(t),
    nFaces_(nF),
    startFace_(sF)
{}
    
boundaryPatchBase::boundaryPatchBase(const word& name, const dictionary& dict)
:
    name_(name),
    type_(),
    nFaces_(),
    startFace_()
{
    word type(dict.lookup("type"));
    type_ = type;
    nFaces_ = readLabel(dict.lookup("nFaces"));
    startFace_ = readLabel(dict.lookup("startFace"));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const boundaryPatchBase& wpb)
{
    wpb.write(os);
    os.check("Ostream& operator<<(Ostream& f, const boundaryPatchBase& wpb");
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
