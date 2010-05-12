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

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamAddToSelection_H
#define vtkPV3FoamAddToSelection_H

// FOAM includes
#include "SortableList.H"

// VTK includes
#include "vtkDataArraySelection.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::vtkPV3Foam::addToSelection
(
    vtkDataArraySelection *arraySelection,
    const IOobjectList& objectLst,
    const string& suffix
)
{
    IOobjectList objects(objectLst.lookupClass(Type::typeName));

    SortableList<word> objectNames(objects.size());

    label count = 0;
    forAllConstIter(IOobjectList, objects, iter)
    {
        objectNames[count++] = iter()->name();
    }

    if (count)
    {
        objectNames.sort();

        forAll (objectNames, objI)
        {
            if (debug)
            {
                Info<<"    addToSelection<" << Type::typeName << "> to GUI "
                    << ":" << objectNames[objI] << endl;
            }

            if (suffix.size())
            {
                arraySelection->AddArray
                (
                    (objectNames[objI] + suffix).c_str()
                );
            }
            else
            {
                arraySelection->AddArray
                (
                    objectNames[objI].c_str()
                );
            }
        }
    }

    return count;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
