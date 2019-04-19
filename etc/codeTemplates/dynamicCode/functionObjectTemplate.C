/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "functionObjectTemplate.H"
#include "functionObject.H"
#include "foamTime.H"
#include "fvCFD.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(${typeName}FunctionObject, 0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const objectRegistry& ${typeName}FunctionObject::obr() const
{
    return obr_;
}


const fvMesh& ${typeName}FunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}FunctionObject::${typeName}FunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool
)
:
    name_(name),
    obr_(obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}FunctionObject::~${typeName}FunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}FunctionObject::read(const dictionary& dict)
{
    if (${verbose:-false})
    {
        Info<<"read ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeRead}
//}}} end code
}


void ${typeName}FunctionObject::execute()
{
    if (${verbose:-false})
    {
        Info<<"execute ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeExecute}
//}}} end code
}


void ${typeName}FunctionObject::end()
{
    if (${verbose:-false})
    {
        Info<<"end ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeEnd}
//}}} end code
}


void ${typeName}FunctionObject::timeSet()
{
    if (${verbose:-false})
    {
        Info<<"timeSet ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin codeTime
    ${codeTimeSet}
//}}} end code
}


void ${typeName}FunctionObject::write()
{
    if (${verbose:-false})
    {
        Info<<"write ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${code}
//}}} end code
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
