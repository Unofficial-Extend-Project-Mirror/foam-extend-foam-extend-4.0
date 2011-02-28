/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany


\*---------------------------------------------------------------------------*/

#include "nasaHeatCapacityPolynomial.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::nasaHeatCapacityPolynomial<equationOfState>::nasaHeatCapacityPolynomial(Istream& is)
:
    equationOfState(is),
    a1_(readScalar(is)),
    a2_(readScalar(is)),
    a3_(readScalar(is)),
    a4_(readScalar(is)),
    a5_(readScalar(is)),
    a6_(readScalar(is)),
    a7_(readScalar(is)),
    //values for some need terms at std
        e0_std(e0(this->Tstd)),
	s0_std(s0(this->Tstd)),
        integral_p_dv_std(this->integral_p_dv(this->rhostd(),this->Tstd)),
        integral_dpdT_dv_std(this->integral_dpdT_dv(this->rhostd(),this->Tstd)),
    // cp @ STD (needed to limit cp for stability
        cp_std(this->cp_nonLimited(this->rhostd(),this->Tstd)) 	
{
    is.check("nasaHeatCapacityPolynomial::nasaHeatCapacityPolynomial(Istream& is)"); 	
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const nasaHeatCapacityPolynomial<equationOfState>& ct
)
{
    os  << static_cast<const equationOfState&>(ct) << tab
        << ct.a1_ << tab<< ct.a2_ << tab << ct.a3_ << tab << ct.a4_ << tab << ct.a5_ << tab << ct.a6_ << tab << ct.a7_ ;
    os.check("Ostream& operator<<(Ostream& os, const nasaHeatCapacityPolynomial& ct)");
    return os;
}


// ************************************************************************* //
