/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Foam::chemistrySolver

Description
    An abstract base class for solving chemistry

SourceFiles
    chemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef chemistrySolver_H
#define chemistrySolver_H

#include "ODEChemistryModel.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CompType, class ThermoType>
class chemistrySolver;

/*---------------------------------------------------------------------------*\
                      Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class chemistrySolver
{
protected:

    // Protected data

        //- Reference to the chemistry model
        ODEChemistryModel<CompType, ThermoType>& model_;

        //- Name of the chemistry solver
        const word name_;


public:

        //- Runtime type information
        TypeName("chemistrySolver");


        // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistrySolver,
            dictionary,
            (
                ODEChemistryModel<CompType, ThermoType>& model,
                const word& modelName
            ),
            (model, modelName)
        );


    // Constructors

        //- Construct from components
        chemistrySolver
        (
            ODEChemistryModel<CompType, ThermoType>& model,
            const word& modelName
        );


    //- Selector
    static autoPtr<chemistrySolver> New
    (
        ODEChemistryModel<CompType, ThermoType>& model,
        const word& compTypeName,
        const word& thermoTypeName
    );


    //- Destructor
    virtual ~chemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual scalar solve
        (
            scalarField &c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistrySolver(Comp, Thermo)                                     \
                                                                              \
    typedef chemistrySolver<Comp, Thermo>                                     \
        chemistrySolver##Comp##Thermo;                                        \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        chemistrySolver##Comp##Thermo,                                        \
        "chemistryModel<"#Comp","#Thermo">",                                  \
        0                                                                     \
    );                                                                        \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        chemistrySolver##Comp##Thermo,                                        \
        dictionary                                                            \
    );


#define makeChemistrySolverType(SS, Comp, Thermo)                             \
                                                                              \
    typedef SS<Comp, Thermo> SS##Comp##Thermo;                                \
                                                                              \
    defineTemplateTypeNameAndDebugWithName                                    \
    (                                                                         \
        SS##Comp##Thermo,                                                     \
        #SS"<"#Comp","#Thermo">",                                             \
        0                                                                     \
    );                                                                        \
                                                                              \
    chemistrySolver<Comp, Thermo>::                                           \
        adddictionaryConstructorToTable<SS<Comp, Thermo> >                    \
            add##SS##Comp##Thermo##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistrySolver.C"
#   include "newChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
