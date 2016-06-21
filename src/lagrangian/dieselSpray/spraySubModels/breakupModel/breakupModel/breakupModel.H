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
    Foam::breakupModel

Description
    Base class for selecting a breakup model

\*---------------------------------------------------------------------------*/

#ifndef breakupModel_H
#define breakupModel_H

#include "IOdictionary.H"
#include "spray.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class breakupModel Declaration
\*---------------------------------------------------------------------------*/

class breakupModel
{

protected:

    // Protected data

        const dictionary& dict_;

        spray& spray_;
        Random& rndGen_;

        Switch includeOscillation_;

        dictionary TABcoeffsDict_;
        scalar y0_;
        scalar yDot0_;
        scalar TABComega_;
        scalar TABCmu_;
        scalar TABWeCrit_;


public:

    //- Runtime type information
        TypeName("breakupModel");

    // Declare runtime constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            breakupModel,
            dictionary,
            (
                const dictionary& dict,
                spray& sm
            ),
            (dict, sm)
        );


    // Constructors

        //- Construct from components
        breakupModel
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        virtual ~breakupModel();


    // Selector

        static autoPtr<breakupModel> New
        (
            const dictionary& dict,
            spray& sm
        );


    // Member Functions

        //- Return the initial deviation from sphericity
        inline scalar y0() const
        {
            return y0_;
        }

        //- Return the initial rate of deviation from sphericity
        inline scalar yDot0() const
        {
            return yDot0_;
        }

        void updateParcelProperties
        (
            parcel& p,
            const scalar deltaT,
            const vector& Ug,
            const liquidMixture& fuels
        ) const;

        virtual void breakupParcel
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

