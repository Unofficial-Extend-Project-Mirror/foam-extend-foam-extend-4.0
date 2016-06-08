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

#include "solidSolver.H"
#include "foamTime.H"
#include "nonLinearGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidSolver> Foam::solidSolver::New(const fvMesh& mesh)
{
    word solidSolverTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the flow is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary solidProperties
        (
            IOobject
            (
                "solidProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        solidProperties.lookup("solidSolver") >> solidSolverTypeName;

        // Add solidMechanics dictionary with nonLinear entry
        if (!mesh.solutionDict().found("solidMechanics"))
        {
            dictionary stressedFoamDict;
            nonLinearGeometry::nonLinearType nonLinear =
                nonLinearGeometry::nonLinearNames_.read
                (
                    solidProperties.subDict(solidSolverTypeName + "Coeffs")
                    .lookup("nonLinear")
                );

            stressedFoamDict.add
                (
                    "nonLinear",
                    nonLinearGeometry::nonLinearNames_[nonLinear]
                );

            const_cast<fvSolution&>(mesh.solutionDict())
           .add("solidMechanics", stressedFoamDict);
        }
        else
        {
            nonLinearGeometry::nonLinearType nonLinear =
                nonLinearGeometry::nonLinearNames_.read
                (
                    solidProperties.subDict(solidSolverTypeName + "Coeffs")
                    .lookup("nonLinear")
                );

            const_cast<dictionary&>
            (
                mesh.solutionDict().subDict("solidMechanics")
            ).set("nonLinear", nonLinearGeometry::nonLinearNames_[nonLinear]);
        }
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solidSolverTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidSolver::New(const fvMesh&)"
        )   << "Unknown solidSolver type " << solidSolverTypeName
            << endl << endl
            << "Valid solidSolver types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidSolver>(cstrIter()(mesh));
}


// ************************************************************************* //
