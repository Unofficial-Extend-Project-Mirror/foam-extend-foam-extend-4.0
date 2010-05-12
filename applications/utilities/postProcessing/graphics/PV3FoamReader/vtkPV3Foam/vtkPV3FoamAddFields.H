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
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamAddFields_H
#define vtkPV3FoamAddFields_H

// FOAM includes
#include "SortableList.H"

// VTK includes
#include "vtkDataArraySelection.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
Foam::label Foam::vtkPV3Foam::addFields
(
    vtkDataArraySelection *fieldSelection,
    const IOobjectList& objects,
    const string& suffix
)
{
    IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

    SortableList<word> fields(fieldObjects.size());

    label count = 0;
    forAllConstIter(IOobjectList, fieldObjects, iter)
    {
        fields[count++] = iter()->name();
    }

    if (count)
    {
        fields.sort();

        forAll(fields, fieldI)
        {
            if (debug)
            {
                Info<<"add region " << GeoField::typeName
                    << ":" << fields[fieldI] << endl;
            }

            if (suffix.size())
            {
                fieldSelection->AddArray
                (
                    (fields[fieldI] + suffix).c_str()
                );
            }
            else
            {
                fieldSelection->AddArray
                (
                    fields[fieldI].c_str()
                );
            }
        }
    }

    return count;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
