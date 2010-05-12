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

Class
    Foam::combustionMixture

Description
    combustion mixture

SourceFiles
    combustionMixture.C

\*---------------------------------------------------------------------------*/

#ifndef combustionMixture_H
#define combustionMixture_H

#include "volFields.H"
#include "PtrList.H"
#include "speciesTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class combustionMixture Declaration
\*---------------------------------------------------------------------------*/

class combustionMixture
{

protected:

    // Protected data

        //- The table of species
        speciesTable species_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


public:

    // Constructors

        //- Construct from dictionary and mesh
        combustionMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&
        );


    // Destructor

        virtual ~combustionMixture()
        {}


    // Member functions

        //- Return the mass-fraction fields
        PtrList<volScalarField>& Y()
        {
            return Y_;
        }

        //- Return the const mass-fraction fields
        const PtrList<volScalarField>& Y() const
        {
            return Y_;
        }

        //- Return the mass-fraction field for a specie given by index
        volScalarField& Y(const label i)
        {
            return Y_[i];
        }

        //- Return the const mass-fraction field for a specie given by index
        const volScalarField& Y(const label i) const
        {
            return Y_[i];
        }

        //- Return the mass-fraction field for a specie given by name
        volScalarField& Y(const word& specieName)
        {
            return Y_[species_[specieName]];
        }

        //- Return the const mass-fraction field for a specie given by name
        const volScalarField& Y(const word& specieName) const
        {
            return Y_[species_[specieName]];
        }

        //- does the mixture include this specie
        bool contains(const word& specieName) const
        {
            return species_.contains(specieName);
        }

        scalar fres(const scalar ft, const scalar stoicRatio) const
        {
            return max(ft - (1.0 - ft)/stoicRatio, 0.0);
        }

        tmp<volScalarField> fres
        (
            const volScalarField& ft,
            const dimensionedScalar& stoicRatio
        ) const
        {
            return max(ft - (scalar(1) - ft)/stoicRatio.value(), scalar(0));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
