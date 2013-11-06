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

#include "fvc.H"
#include "leastSquaresGrad.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Store gradients of fields on the mesh prior to topology changes
template <class Type, class gradType>
void topoMapper::storeGradients
(
    HashTable<autoPtr<gradType> >& gradTable
) const
{
    // Define a few typedefs for convenience
    typedef typename outerProduct<vector, Type>::type gCmptType;
    typedef GeometricField<Type, fvPatchField, volMesh> volType;
    typedef GeometricField<gCmptType, fvPatchField, volMesh> gVolType;

    // Fetch all fields from registry
    HashTable<const volType*> fields
    (
        mesh_.objectRegistry::lookupClass<volType>()
    );

    // Store old-times before gradient computation
    forAllIter(typename HashTable<const volType*>, fields, fIter)
    {
        fIter()->storeOldTimes();
    }

    forAllConstIter(typename HashTable<const volType*>, fields, fIter)
    {
        const volType& field = *fIter();

        // Compute the gradient.
        tmp<gVolType> tGrad;

        // If the fvSolution dictionary contains an entry,
        // use that, otherwise, default to leastSquares
        word gradName("grad(" + field.name() + ')');

        if (mesh_.schemesDict().subDict("gradSchemes").found(gradName))
        {
            tGrad = fvc::grad(field);
        }
        else
        {
            tGrad = fv::leastSquaresGrad<Type>(mesh_).grad(field);
        }

        // Make a new entry, but don't register the field.
        gradTable.insert
        (
            field.name(),
            autoPtr<gradType>
            (
                new gradType
                (
                    IOobject
                    (
                        tGrad().name(),
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    tGrad()
                )
            )
        );
    }
}


} // End namespace Foam

// ************************************************************************* //
