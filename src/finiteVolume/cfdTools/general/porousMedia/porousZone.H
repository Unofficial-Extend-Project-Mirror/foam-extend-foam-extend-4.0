/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::porousZone

Description
    Porous zone definition based on cell zones.

    Porous zone definition based on cell zones and parameters obtained from a
    control dictionary constructed from the given stream. The orientation of
    the porous region is defined with the same notation as a coordinateSystem,
    but only a Cartesian coordinate system is valid.

    Implemented porosity models:

    powerLaw (@e C0 and @e C1 parameters)
    @f[
        S = - \rho C_0 |U|^{(C_1 - 1)/2} U
    @f]

    In the porous zone a heat source using the dual stream heat exchanger
    approach:

    heatTransfer
    @f[
        Maux = Mass flow of aux fluid                         // [kg/s]
        Qaux = Heat to be eliminated                          // [W]
        Caux = Specific heat capacity of aux fluid            // [J / kg K]
        Taux = initial guess for aux fluid inlet temperature  // [K]
    @f]

    Darcy-Forchheimer (@e d and @e f parameters)
    @f[
        S = - (\mu \, d + \frac{\rho |U|}{2} \, f) U
    @f]

    Since negative Darcy/Forchheimer parameters are invalid, they can be used
    to specify a multiplier (of the max component).

    The porousZones method porousZones::ddt() mirrors the normal fvm::ddt()
    method, but accounts for the effective volume of the cells.

See Also
    porousZones and coordinateSystems

SourceFiles
    porousZone.C
    porousZoneTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef porousZone_H
#define porousZone_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"

#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

// Forward declaration of friend functions and operators
class porousZone;

Ostream& operator<<(Ostream&, const porousZone&);


/*---------------------------------------------------------------------------*\
                        Class porousZone Declaration
\*---------------------------------------------------------------------------*/

class porousZone
{
    // Private data

        //- Name of this zone
        word name_;

        //- Reference to the finite volume mesh this zone is part of
        const fvMesh& mesh_;

        //- Dictionary containing the parameters
        dictionary dict_;

        //- Cell zone ID
        label cellZoneID_;

        //- Coordinate system used for the zone (Cartesian)
        coordinateSystem coordSys_;

        //- Porosity of the zone (0 < porosity <= 1)
        //  Placeholder for treatment of temporal terms.
        //  Currently unused.
        scalar porosity_;

        //- Power Law coefficient C0
        scalar C0_;

        //- Power Law coefficient C1
        scalar C1_;

        //- Darcy coefficient
        dimensionedTensor D_;

        //- Forchheimer coefficient
        dimensionedTensor F_;

        //- Density of primary fluid
        scalar rhoPri_;

        //- Specific heat capacity of primary fluid
        scalar Cpri_;
    
        //- Mass flow of aux fluid
        scalar Maux_;

        //- Heat to be eliminated
        scalar Qaux_;

        //- Specific heat capacity of coolant
        scalar Caux_;

        //- Inlet temperature of aux fluid
        mutable scalar Taux_;

        //- Epsilon for Single Effectiveness HX
        scalar Qepsilon_;

        //- Relaxation factor for Taux
        scalar TauxRelax_;


        // Geometrical description of the porous heat exchanger

            //- Number of inlet cells for aux fluid
            label nCellsAuxInlet_;

            //- First inlet cell for aux fluid
            label firstCell_;

            //- Inlet direction for auxiliary fluid
            vector auxUnitVector_;

            //- Number of vertical cells for auxiliary fluid
            label nVerticalCells_;


    // Private Member Functions

        //- Adjust negative resistance values to be multiplier of max value
        static void adjustNegativeResistance(dimensionedVector& resist);

        //- Power-law resistance
        template<class RhoFieldType>
        void addPowerLawResistance
        (
            scalarField& Udiag,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            scalarField& Udiag,
            vectorField& Usource,
            const labelList& cells,
            const scalarField& V,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Power-law resistance
        template<class RhoFieldType>
        void addPowerLawResistance
        (
            tensorField& AU,
            const labelList& cells,
            const RhoFieldType& rho,
            const vectorField& U
        ) const;

        //- Viscous and inertial resistance
        template<class RhoFieldType>
        void addViscousInertialResistance
        (
            tensorField& AU,
            const labelList& cells,
            const RhoFieldType& rho,
            const scalarField& mu,
            const vectorField& U
        ) const;

        //- Heat source
        template<class RhoFieldType>
        void addHeatSource
        (
            const scalarField& Macro,
            const scalarField& posFlux,
            scalarField& QSource,
            scalarField& Qaux,
            const scalarField& T,
            scalarField& Taux,
            const vectorField& U,
            const RhoFieldType& rho
        ) const;


        //- Disallow default bitwise copy construct
        porousZone(const porousZone&);

        //- Disallow default bitwise assignment
        void operator=(const porousZone&);


public:

    // Constructors

        //- Construct from components
        porousZone(const word& name, const fvMesh&, const dictionary&);

        //- Return clone
        autoPtr<porousZone> clone() const
        {
            notImplemented("autoPtr<porousZone> clone() const");
            return autoPtr<porousZone>(nullptr);
        }

        //- Return pointer to new porousZone created on freestore from Istream
        class iNew
        {
            //- Reference to the finite volume mesh this zone is part of
            const fvMesh& mesh_;

        public:

            iNew(const fvMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<porousZone> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);

                return autoPtr<porousZone>(new porousZone(name, mesh_, dict));
            }
        };


    //- Destructor
    virtual ~porousZone()
    {}


    // Member Functions

        // Access

            //- cellZone name
            const word& zoneName() const
            {
                return name_;
            }

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return cellZone number
            label zoneId() const
            {
                return cellZoneID_;
            }

            //- Return cellZone
            const cellZone& zone() const
            {
                return mesh_.cellZones()[cellZoneID_];
            }

            //- Dictionary values used for the porousZone
            const dictionary& dict() const
            {
                return dict_;
            }

            //- Return coordinate system
            const coordinateSystem& coordSys() const
            {
                return coordSys_;
            }

            //- Return origin
            const point& origin() const
            {
                return coordSys_.origin();
            }

            //- Return axis
            vector axis() const
            {
                return coordSys_.axis();
            }

            //- Return porosity
            scalar porosity() const
            {
                return porosity_;
            }

            //- Edit access to porosity
            scalar& porosity()
            {
                return porosity_;
            }


        //- Modify time derivative elements according to porosity
        template<class Type>
        void modifyDdt(fvMatrix<Type>&) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the momentum equation
        void addResistance(fvVectorMatrix& UEqn) const;

        //- Add the viscous and inertial resistance force contribution
        //  to the tensorial diagonal.
        //  Optionally correct the processor BCs of AU.
        void addResistance
        (
            const fvVectorMatrix& UEqn,
            volTensorField& AU,
            bool correctAUprocBC = true
        ) const;

        //- Add the heat source in enthalpy equation
        void addHeatResistance
        (
            fvScalarMatrix& hTEqn,
            const volScalarField& T,
            volScalarField& Taux,
            volScalarField& Qaux,
            const volVectorField& U,
            const volScalarField& Macro,
            const volScalarField& posFlux
        ) const;

        //- Order cells for Dual Stream model
        void macroCellOrder
        (
            volScalarField& Taux,
            volScalarField& Macro,
            volScalarField& posFlux,
            const surfaceScalarField& phi
        ) const;

        //- Write the porousZone dictionary
        virtual void writeDict(Ostream&, bool subDict = true) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const porousZone&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "porousZoneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
