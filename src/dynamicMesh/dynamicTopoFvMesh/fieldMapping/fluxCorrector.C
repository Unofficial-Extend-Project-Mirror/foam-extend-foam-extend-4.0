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

Class
    fluxCorrector

Description
    Implementation of the fluxCorrector base class

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "fluxCorrector.H"
#include "dlLibraryTable.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(fluxCorrector, 0);
defineRunTimeSelectionTable(fluxCorrector, mesh);

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<fluxCorrector> fluxCorrector::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    // Check if an optional entry was specified
    if (dict.found("fluxCorrector"))
    {
        word correctorTypeName(dict.lookup("fluxCorrector"));

        // Open any supplied libraries in dictionary
        dlLibraryTable::open
        (
            dict,
            "fluxCorrectorLibs",
            meshConstructorTablePtr_
        );

        if (!meshConstructorTablePtr_)
        {
            FatalErrorIn
            (
                "autoPtr<fluxCorrector> fluxCorrector::New"
                "(const fvMesh& mesh, const dictionary& dict)"
            )   << "fluxCorrector table is empty"
                << exit(FatalError);
        }

        correctorTypeName = word(dict.lookup("fluxCorrector"));

        meshConstructorTable::iterator cstrIter =
            meshConstructorTablePtr_->find(correctorTypeName);

        if (cstrIter == meshConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "autoPtr<fluxCorrector> fluxCorrector::New"
                "(const fvMesh& mesh, const dictionary& dict)"
            )   << "Unknown fluxCorrector type " << correctorTypeName
                << endl << endl
                << "Valid fluxCorrector types are: " << endl
                << meshConstructorTablePtr_->toc()
                << exit(FatalError);
        }

        return autoPtr<fluxCorrector>(cstrIter()(mesh, dict));
    }

    // Return the default fluxCorrector
    return autoPtr<fluxCorrector>(new fluxCorrector(mesh, dict));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return reference to mesh
const fvMesh& fluxCorrector::mesh() const
{
    return mesh_;
}


//- Return reference to dictionary
const dictionary& fluxCorrector::dict() const
{
    return dict_;
}


//- Is flux-correction required?
bool fluxCorrector::required() const
{
    // Support for cases which do not involve a flow-solver.
    // Notify the user, just in case.
    Info << " ~~~ No flux correction ~~~ " << endl;

    return false;
}


//- Interpolate fluxes to a specified list of faces
void fluxCorrector::interpolateFluxes(const labelList& faces) const
{
    if (required())
    {
        // Throw an error stating that the scheme hasn't been implemented
        notImplemented
        (
            "void fluxCorrector::interpolateFluxes"
            "(const labelList& faces) const"
        );
    }
}


//- Update fluxes in the registry, if required
void fluxCorrector::updateFluxes() const
{
    if (required())
    {
        // Throw an error stating that the scheme hasn't been implemented
        notImplemented("void fluxCorrector::updateFluxes() const");
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
void fluxCorrector::operator=(const fluxCorrector& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("fluxCorrector::operator=(const fluxCorrector&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

} // End namespace Foam

// ************************************************************************* //
