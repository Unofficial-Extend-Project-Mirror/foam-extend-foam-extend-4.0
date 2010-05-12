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
    Foam::CH3OH

Description
    methanol

SourceFiles
    CH3OH.C

\*---------------------------------------------------------------------------*/

#ifndef CH3OH_H
#define CH3OH_H

#include "liquid.H"
#include "NSRDSfunc0.H"
#include "NSRDSfunc1.H"
#include "NSRDSfunc2.H"
#include "NSRDSfunc3.H"
#include "NSRDSfunc4.H"
#include "NSRDSfunc5.H"
#include "NSRDSfunc6.H"
#include "NSRDSfunc7.H"
#include "NSRDSfunc14.H"
#include "APIdiffCoefFunc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CH3OH Declaration
\*---------------------------------------------------------------------------*/

class CH3OH
:
    public liquid
{
    // Private data

        NSRDSfunc5 rho_;
        NSRDSfunc1 pv_;
        NSRDSfunc6 hl_;
        NSRDSfunc0 cp_;
        NSRDSfunc0 h_;
        NSRDSfunc7 cpg_;
        NSRDSfunc4 B_;
        NSRDSfunc1 mu_;
        NSRDSfunc2 mug_;
        NSRDSfunc0 K_;
        NSRDSfunc2 Kg_;
        NSRDSfunc6 sigma_;
        APIdiffCoefFunc D_;


public:

    //- Runtime type information
    TypeName("CH3OH");


    // Constructors

        //- Construct null
        CH3OH()
        :
            liquid(32.042, 512.58, 8.0959e+6, 0.1178, 0.224, 175.47, 1.054e-1, 337.85, 5.6706e-30, 0.5656, 2.9523e+4),
            rho_(73.952936, 0.27192, 512.58, 0.2331),
            pv_(109.93, -7471.3, -13.988, 0.015281, 1),
            hl_(512.58, 1644716.30984333, 0.3766, 0, 0, 0),
            cp_(3358.09250358904, -11.8781599151114, 0.0305536483365583, 0, 0, 0),
            // NN: enthalpy, h_, is not used in the sprayModel.
            // For consistency, the enthalpy is derived from hlat and hl.
            // It is, however, convenient to have it available.
            h_(-8190474.32066862, 3358.09250358904, -5.93907995755571, 0.0101845494455194, 0, 0),
            cpg_(1226.9521253355, 2772.92303851195, 1963, 1733.66206853505, 909.6),
            B_(-0.0199737844079645, 19.3496036452157, -3342487.98452032, 2.40808938268523e+19, -6.85787404032208e+21),
            mu_(-7.288, 1065.3, -0.6657, 0, 0),
            mug_(3.0663e-07, 0.69655, 205, 0),
            K_(0.2837, -0.000281, 0, 0, 0, 0),
            Kg_(-7.763, 1.0279, -74360000.0, 6770000000.0),
            sigma_(512.58, 0.056, -0.00014583, 1.08e-07, 0, 0),
            D_(147.18, 20.1, 32.042, 28) // NN: Same as nHeptane
        {}
        CH3OH
        (
            const liquid& l,
            const NSRDSfunc5& density,
            const NSRDSfunc1& vapourPressure,
            const NSRDSfunc6& heatOfVapourisation,
            const NSRDSfunc0& heatCapacity,
            const NSRDSfunc0& enthalpy,
            const NSRDSfunc7& idealGasHeatCapacity,
            const NSRDSfunc4& secondVirialCoeff,
            const NSRDSfunc1& dynamicViscosity,
            const NSRDSfunc2& vapourDynamicViscosity,
            const NSRDSfunc0& thermalConductivity,
            const NSRDSfunc2& vapourThermalConductivity,
            const NSRDSfunc6& surfaceTension,
            const APIdiffCoefFunc& vapourDiffussivity
        )
        :
            liquid(l),
            rho_(density),
            pv_(vapourPressure),
            hl_(heatOfVapourisation),
            cp_(heatCapacity),
            h_(enthalpy),
            cpg_(idealGasHeatCapacity),
            B_(secondVirialCoeff),
            mu_(dynamicViscosity),
            mug_(vapourDynamicViscosity),
            K_(thermalConductivity),
            Kg_(vapourThermalConductivity),
            sigma_(surfaceTension),
            D_(vapourDiffussivity)
        {}

        //- Construct from Istream
        CH3OH(Istream& is)
        :
            liquid(is),
            rho_(is),
            pv_(is),
            hl_(is),
            cp_(is),
            h_(is),
            cpg_(is),
            B_(is),
            mu_(is),
            mug_(is),
            K_(is),
            Kg_(is),
            sigma_(is),
            D_(is)
        {}


    // Member Functions

        //- Liquid density [kg/m^3]
        scalar rho(scalar p, scalar T) const
        {
            return rho_.f(p, T);
        }

        //- Vapour pressure [Pa]
        scalar pv(scalar p, scalar T) const
        {
            return pv_.f(p, T);
        }

        //- Heat of vapourisation [J/kg]
        scalar hl(scalar p, scalar T) const
        {
            return hl_.f(p, T);
        }

        //- Liquid heat capacity [J/(kg K)]
        scalar cp(scalar p, scalar T) const
        {
            return cp_.f(p, T);
        }

        //- Liquid Enthalpy [J/(kg)]
        scalar h(scalar p, scalar T) const
        {
            return h_.f(p, T);
        }

        //- Ideal gas heat capacity [J/(kg K)]
        scalar cpg(scalar p, scalar T) const
        {
            return cpg_.f(p, T);
        }

        //- Second Virial Coefficient [m^3/kg]
        scalar B(scalar p, scalar T) const
        {
            return B_.f(p, T);
        }

        //- Liquid viscosity [Pa s]
        scalar mu(scalar p, scalar T) const
        {
            return mu_.f(p, T);
        }

        //- Vapour viscosity [Pa s]
        scalar mug(scalar p, scalar T) const
        {
            return mug_.f(p, T);
        }

        //- Liquid thermal conductivity  [W/(m K)]
        scalar K(scalar p, scalar T) const
        {
            return K_.f(p, T);
        }

        //- Vapour thermal conductivity  [W/(m K)]
        scalar Kg(scalar p, scalar T) const
        {
            return Kg_.f(p, T);
        }

        //- Surface tension [N/m]
        scalar sigma(scalar p, scalar T) const
        {
            return sigma_.f(p, T);
        }

        //- Vapour diffussivity [m2/s]
        scalar D(scalar p, scalar T) const
        {
            return D_.f(p, T);
        }


        //- Write the function coefficients
        void writeData(Ostream& os) const
        {
            liquid::writeData(os); os << nl;
            rho_.writeData(os); os << nl;
            pv_.writeData(os); os << nl;
            hl_.writeData(os); os << nl;
            cp_.writeData(os); os << nl;
            cpg_.writeData(os); os << nl;
            B_.writeData(os); os << nl;
            mu_.writeData(os); os << nl;
            mug_.writeData(os); os << nl;
            K_.writeData(os); os << nl;
            Kg_.writeData(os); os << nl;
            sigma_.writeData(os); os << nl;
            D_.writeData(os); os << endl;
        }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const CH3OH& l)
        {
            l.writeData(os);
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
