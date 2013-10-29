/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidInterface

Description

\*---------------------------------------------------------------------------*/

#include "solidInterface.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

autoPtr<solidInterface> solidInterface::New
(
 const word& name,
 const fvMesh& mesh,
 const constitutiveModel& rheology
 )
{
  Info<< "Selecting solidInterface procedure " << name << endl;

  dictionaryConstructorTable::iterator cstrIter =
    dictionaryConstructorTablePtr_->find(name);
  
    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "solidInterface::New(\n"
            "    const word& name,\n"
            "    const fvMesh& mesh,\n"
            "    const constitutiveModel& rheology" //,\n"
            ")",
	    rheology
        )   << "Unknown solidInterfaceMethod "
            << name << endl << endl
            << "Valid  solidInterfaces are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<solidInterface>(cstrIter()(name, mesh, rheology));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
