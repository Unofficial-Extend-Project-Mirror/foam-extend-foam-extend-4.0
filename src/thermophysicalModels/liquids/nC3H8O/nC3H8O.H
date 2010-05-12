/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   nC3H8O
   \\  /           | Family: nC3H8O
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2000 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    nC3H8O

DESCRIPTION
    propanol 

*/
// ------------------------------------------------------------------------- //

#ifndef nC3H8O_H
#define nC3H8O_H

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class nC3H8O Declaration
\*---------------------------------------------------------------------------*/

class nC3H8O
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
        NSRDSfunc0 sigma_;
        NSRDSfunc1 D_;

public:

    //- Runtime type information
    TypeName("nC3H8O");


    // Constructors

        //- Construct null
        nC3H8O()
        :
            liquid(60.096, 536.71, 5.1696e+6, 0.21853, 0.253, 146.95, 6.5112e-7, 370.35, 5.6039e-30, 0.6279, 2.4557e+4),
            rho_(75.300288, 0.272, 536.71, 0.2494),
            pv_(77.46, -7960, -7.5235, 3e-07, 2),
            hl_(536.71, 1098242.8115016, 0.647, -0.783, 0.613, 0),
            cp_(216.320553780618, 18.5203674121406, -0.0751797124600639, 0.000126464323748669, 0, 0),
            // NN: enthalpy, h_, is not used in the sprayModel.
            // For consistency, the enthalpy is derived from hlat and hl.
            // It is, however, convenient to have it available.
            h_(-5533091.96851587, 216.320553780618, 9.26018370607029, -0.0250599041533546, 3.16160809371672e-05, 0),
            cpg_(961.794462193823, 3467.78487752929, 1542, 2046.72523961661, 649),
            B_(0.000933506389776358, -1.09325079872204, -531649.361022364, -2.32627795527157e+17, -3.81888977635783e+20),
            mu_(0.571, 1521, -2.0894, 0, 0),
            mug_(7.942e-07, 0.5491, 415.8, 0),
            K_(0.204, -0.000169, 0, 0, 0, 0),
            Kg_(-613.84, 0.7927, -1157400000, 0),
            sigma_(0.04533, -6.88e-05, -1.6e-08, 0, 0, 0),
            D_(4.75e-10, 1.75, 0.0, 0.0, 0.0) // NN. same as iC3H8O
        {}
        nC3H8O
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
            const NSRDSfunc0& surfaceTension,
            const NSRDSfunc1& vapourDiffussivity
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
        nC3H8O(Istream& is)
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

        friend Ostream& operator<<(Ostream& os, const nC3H8O& l)
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

