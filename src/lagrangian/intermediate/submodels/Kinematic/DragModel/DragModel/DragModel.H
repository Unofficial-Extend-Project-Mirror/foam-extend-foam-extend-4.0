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
    Foam::DragModel


Description
    Templated drag model class

SourceFiles
    DragModel.C
    NewDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef DragModel_H
#define DragModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class DragModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DragModel
{

    // Private data

        //- The cloud dictionary
        const dictionary& dict_;

        // reference to the owner cloud class
        CloudType& owner_;


public:

    //- Runtime type information
    TypeName("DragModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DragModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct from components
        DragModel
        (
            const dictionary& dict,
            CloudType& owner
        );


    // Destructor

        virtual ~DragModel();


    // Selector

        static autoPtr<DragModel<CloudType> > New
        (
            const dictionary& dict,
            CloudType& cloud
        );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the dictionary
        const dictionary& dict() const;


    // Member Functions

        //- Flag to indicate whether model activates drag model
        virtual bool active() const = 0;

        //- Return drag coefficient
        virtual scalar Cd(const scalar Re) const = 0;

        //- Return linearised coefficient for velocity equation
        //  Drag force per unit particle mass = Cu(U - Up)
        scalar Cu
        (
            const vector& Ur,
            const scalar d,
            const scalar rhoc,
            const scalar rhop,
            const scalar mu
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDragModel(CloudType)                                              \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(DragModel<CloudType>, 0);             \
                                                                              \
    defineTemplateRunTimeSelectionTable(DragModel<CloudType>, dictionary);


#define makeDragModelType(SS, CloudType, ParcelType)                          \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<CloudType<ParcelType> >, 0);       \
                                                                              \
    DragModel<CloudType<ParcelType> >::                                       \
        adddictionaryConstructorToTable<SS<CloudType<ParcelType> > >          \
            add##SS##CloudType##ParcelType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DragModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
