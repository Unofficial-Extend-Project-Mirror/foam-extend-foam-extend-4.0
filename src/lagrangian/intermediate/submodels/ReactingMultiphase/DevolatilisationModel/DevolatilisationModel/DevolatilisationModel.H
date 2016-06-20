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
    Foam::DevolatilisationModel

Description
    Templated devolatilisation model class

SourceFiles
    DevolatilisationModel.C
    NewDevolatilisationModel.C

\*---------------------------------------------------------------------------*/

#ifndef DevolatilisationModel_H
#define DevolatilisationModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class DevolatilisationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DevolatilisationModel
{
protected:

    // Protected data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficient dictionary
        const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("DevolatilisationModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DevolatilisationModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        DevolatilisationModel(CloudType& owner);

        //- Construct from dictionary
        DevolatilisationModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~DevolatilisationModel();


    //- Selector
    static autoPtr<DevolatilisationModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the cloud dictionary
        const dictionary& dict() const;

        //- Return the coefficient dictionary
        const dictionary& coeffDict() const;


    // Member Functions

        //- Flag to indicate whether model activates devolatilisation model
        virtual bool active() const = 0;

        //- Update model
        virtual scalar calculate
        (
            const scalar dt,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalar YVolatile0,
            const scalar YVolatile,
            bool& canCombust
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDevolatilisationModel(CloudType)                                  \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(DevolatilisationModel<CloudType>, 0); \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        DevolatilisationModel<CloudType>,                                     \
        dictionary                                                            \
    );


#define makeDevolatilisationModelThermoType(SS, CloudType, ParcelType, ThermoType)\
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    DevolatilisationModel<CloudType<ParcelType<ThermoType> > >::              \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DevolatilisationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
