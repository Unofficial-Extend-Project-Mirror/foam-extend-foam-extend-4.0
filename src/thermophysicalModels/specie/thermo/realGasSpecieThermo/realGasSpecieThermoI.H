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

#include "realGasSpecieThermo.H"

template<class thermo>
inline Foam::realGasSpecieThermo<thermo>::realGasSpecieThermo
(
    const thermo& sp
)
:
    thermo(sp)
{}



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// using a stabilizing newton solver 
template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::T
(
    scalar f,
    scalar &T0,
    scalar p,
    scalar &rho0,
    scalar (realGasSpecieThermo<thermo>::*F)(const scalar,const scalar) const,
    scalar (realGasSpecieThermo<thermo>::*dFdT)(const scalar,const scalar) const
) const
{

   
    scalar Test = T0;
    scalar Tnew = T0;
    scalar rho=rho0;
    scalar Ttol = T0*tol_;
    int    iter = 0;
int i;

    do
    {
        Test = Tnew;
        
	rho=this->rho(p,Test,rho);


   i=0;
;
      	do
   	{
        Tnew = Test - ((this->*F)(rho,Test) - f)/(this->*dFdT)(rho,Test)/(pow(2,i)); 
	i++;				// if the solve is diverging, the time step is reduced until the solver converges	
         }while((i<20)&&((mag((this->*F)(rho,Tnew) - f) > mag((this->*F)(rho,Test) - f))));	

        if (iter++ > maxIter_)
        {
            FatalErrorIn
            (
                "realGasSpecieThermo<thermo>::T(scalar f, scalar T0, "
                "scalar (realGasSpecieThermo<thermo>::*F)(const scalar) const, "
                "scalar (realGasSpecieThermo<thermo>::*dFdT)(const scalar) const"
                ") const"
            )   << "Maximum number of iterations exceeded"
                << abort(FatalError);
        }


    } while (mag(mag(Tnew) - mag(Test)) > Ttol);		

  rho0=rho;
  T0=Tnew;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::realGasSpecieThermo<thermo>::realGasSpecieThermo
(
    const word& name,
    const realGasSpecieThermo& st
)
:
    thermo(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::g(const scalar rho,  const scalar T ) const
{
    return this->h(rho, this->pReturn(rho,T)) - T*this->s(rho, this->pReturn(rho,T));
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::a(const scalar rho, const scalar T ) const
{
    return this->e(rho,this->pReturn(rho,T)) - T*this->s(rho,this->pReturn(rho,T));
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::Cp( const scalar rho, const scalar T) const
{
    return this->cp(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::Cv( const scalar rho, const scalar T) const
{
    return this->cv(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::H(const scalar rho, const scalar T) const
{		
    return this->h(rho,  T)/this->W();
}





template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::S(const scalar rho, const scalar T) const
{
    return this->s(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::E(const scalar rho, const scalar T) const
{
    return this->e(rho, T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::G(const scalar rho,  const scalar T) const
{
    return this->g(rho,  T)/this->W();
}


template<class thermo>
inline Foam::scalar Foam::realGasSpecieThermo<thermo>::A(const scalar rho,  const scalar T) const
{
    return this->a(rho, T)/this->W();
}




template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::TH
(
    const scalar h,
    scalar &T0,
    const scalar p,
    scalar &rho0
) const
{
   T(h, T0, p,rho0, &realGasSpecieThermo<thermo>::H, &realGasSpecieThermo<thermo>::Cp);
}


template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::TE
(
    const scalar e,
    scalar &T0,
    const scalar p,
    scalar &rho0
)  const
{
    T(e, T0, p,rho0, &realGasSpecieThermo<thermo>::E, &realGasSpecieThermo<thermo>::Cv);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::operator+=
(
    const realGasSpecieThermo<thermo>& st
)
{
    thermo::operator+=(st);
}

template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::operator-=
(
    const realGasSpecieThermo<thermo>& st
)
{
    thermo::operator-=(st);
}

template<class thermo>
inline void Foam::realGasSpecieThermo<thermo>::operator*=(const scalar s)
{
    thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator+
(
    const realGasSpecieThermo<thermo>& st1,
    const realGasSpecieThermo<thermo>& st2
)
{
    return realGasSpecieThermo<thermo>
    (
        static_cast<const thermo&>(st1) + static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator-
(
    const realGasSpecieThermo<thermo>& st1,
    const realGasSpecieThermo<thermo>& st2
)
{
    return realGasSpecieThermo<thermo>
    (
        static_cast<const thermo&>(st1) - static_cast<const thermo&>(st2)
    );
}


template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator*
(
    const scalar s,
    const realGasSpecieThermo<thermo>& st
)
{
    return realGasSpecieThermo<thermo>
    (
        s*static_cast<const thermo&>(st)
    );
}


template<class thermo>
inline Foam::realGasSpecieThermo<thermo> Foam::operator==
(
    const realGasSpecieThermo<thermo>& st1,
    const realGasSpecieThermo<thermo>& st2
)
{
    return st2 - st1;
}


// ************************************************************************* //
