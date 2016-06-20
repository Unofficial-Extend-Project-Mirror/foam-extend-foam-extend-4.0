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
    Foam::atomizationModel

Description
    Base class for selecting the atomization model

\*---------------------------------------------------------------------------*/

#ifndef atomizationModel_H
#define atomizationModel_H

#include "IOdictionary.H"
#include "spray.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class atomizationModel Declaration
\*---------------------------------------------------------------------------*/

class atomizationModel
{

protected:

    // Protected data

        const dictionary& dict_;
        spray& spray_;
        Random& rndGen_;

public:

    //- Runtime type information
        TypeName("atomizationModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            atomizationModel,
            dictionary,
            (
                const dictionary& dict,
                spray& sm
            ),
            (dict, sm)
        );


    // Constructors

        //- Construct from components
        atomizationModel
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        virtual ~atomizationModel();


    // Selector

        static autoPtr<atomizationModel> New
        (
            const dictionary& dict,
            spray& sm
        );


    // Member Functions

        virtual void atomizeParcel
        (
            parcel& parcel,
            const scalar deltaT,
            const vector& vel,
            const liquidMixture& fuels
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

