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
    Foam::ETAB

Description
    The Enhanced %TAB model.

    Described in the papers below.
    @verbatim
    F.X. Tanner
        "Liquid Jet Atomization and Droplet Breakup Modeling of
        Non-Evaporating Diesel Fuel Sprays"
        SAE 970050,
        SAE Transactions: Journal of Engines, Vol 106, Sec 3 pp 127-140

    F.X. Tanner and G. Weisser
        "Simulation of Liquid Jet Atomization for
        Fuel Sprays by Means of Cascade Drop Breakup Model"
        SAE 980808
        SAE Technical Paper Series
    @endverbatim

See Also
    The TAB model

SourceFiles
    ETAB.C

\*---------------------------------------------------------------------------*/

#ifndef ETAB_H
#define ETAB_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ETAB Declaration
\*---------------------------------------------------------------------------*/

class ETAB
:
    public breakupModel
{
private:

    // Private data

        dictionary coeffsDict_;

        // model constants

        // Cmu_ and Comega_ are the same as in the TAB model
        scalar Cmu_;
        scalar Comega_;

        scalar k1_;
        scalar k2_;
        scalar WeCrit_;
        scalar WeTransition_;
        scalar AWe_;

public:


    //- Runtime type information
    TypeName("ETAB");


    // Constructors

        //- Construct from components
        ETAB
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~ETAB();


    // Member Functions

        void breakupParcel
        (
            parcel& parcel,
            const scalar deltaT,
            const vector& Ug,
            const liquidMixture& fuels
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
