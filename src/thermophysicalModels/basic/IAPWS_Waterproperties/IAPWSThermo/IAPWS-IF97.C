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

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#include "IAPWS-IF97.H"

//CL: calculated all (minimal) needed properties for a given pressure and enthalpy
void Foam::calculateProperties_ph
(
    scalar &p,
    scalar &h,
    scalar &T,
    scalar &rho,
    scalar &psi,
    scalar &drhodh,
    scalar &mu,
    scalar &alpha
)
{
    SteamState S;

    // CL: vapor mass fraction is also calculated in calculateProperties_h
    // CL: in this fuction, x is a dummy variable and x is not return to IAPWSThermo.C
    scalar x;

    S=freesteam_set_ph(p,h);
    calculateProperties_h(S,p,h,T,rho,psi,drhodh,mu,alpha,x);
}


//CL: calculated all (minimal) needed properties + the vapor mass fraction for a given pressure and enthalpy
void Foam::calculateProperties_ph
(
    scalar &p,
    scalar &h,
    scalar &T,
    scalar &rho,
    scalar &psi,
    scalar &drhodh,
    scalar &mu,
    scalar &alpha,
    scalar &x
)
{
    SteamState S;

    S=freesteam_set_ph(p,h);
    calculateProperties_h(S,p,h,T,rho,psi,drhodh,mu,alpha,x);
}


//CL: calculated all (minimal) needed properties for a given pressure and temperature
void Foam::calculateProperties_pT
(
    scalar &p,
    scalar &T,
    scalar &h,
    scalar &rho,
    scalar &psi,
    scalar &drhodh,
    scalar &mu,
    scalar &alpha
)
{
    SteamState S;

    // CL: vapor mass fraction is also calculated in calculateProperties_h
    // CL: in this fuction, x is a dummy variable and x is not return to IAPWSThermo.
    scalar x;

    S=freesteam_set_pT(p,T);
    calculateProperties_h(S,p,h,T,rho,psi,drhodh,mu,alpha,x);
}


//CL: calculated all (minimal) needed properties + the vapor mass fraction for a given pressure and temperature
void Foam::calculateProperties_pT
(
    scalar &p,
    scalar &T,
    scalar &h,
    scalar &rho,
    scalar &psi,
    scalar &drhodh,
    scalar &mu,
    scalar &alpha,
    scalar &x
)
{
    SteamState S;

    S=freesteam_set_pT(p,T);
    calculateProperties_h(S,p,h,T,rho,psi,drhodh,mu,alpha,x);
}


//CL: calculated the properties --> this function is called by the functions above
//CL: does not calulated the internal energy, if this is needed e.g. for sonicFoam
//CL: the function has to be changed a little bit
void Foam::calculateProperties_h
(
    SteamState S,
    scalar &p,
    scalar &h,
    scalar &T,
    scalar &rho,
    scalar &psi,
    scalar &drhodh,
    scalar &mu,
    scalar &alpha,
    scalar &x
)
{
    label region;
    scalar kappa,lambda,cp,beta;

    region=freesteam_region(S);

    //CL:Liquid phase
    if (region==1)
    {
        p=S.R1.p;
        T=S.R1.T;
        rho=1/freesteam_region1_v_pT(S.R1.p,S.R1.T);
        h=freesteam_region1_h_pT(S.R1.p,S.R1.T);
        x=0;

        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        //Cl: note: in FreeStream, kappa=1/V*(dV/dP)_T=const is called kappaT (in this region)
        kappa=freesteam_region1_kappaT_pT(S.R1.p,S.R1.T);
        beta=freesteam_region1_alphav_pT(S.R1.p,S.R1.T);
        cp=freesteam_region1_cp_pT(S.R1.p,S.R1.T);

        //CL: getting derivatives using Bridgmans table
        //CL: psi=(drho/dp)_h=const
        //CL: drhodh=(drho/dh)_p=const
        psi=-((T*beta*beta-beta)/cp-kappa*rho);
        drhodh=-rho*beta/cp;

        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
        lambda=freesteam_k_rhoT(rho,T);
        alpha=lambda/cp; //Cl: Important info -->alpha= thermal diffusivity time density
    }
    //CL:vapor phase
    else if (region==2)
    {
        p=S.R2.p;
        T=S.R2.T;
        rho=1/freesteam_region2_v_pT(S.R2.p,S.R2.T);
        h=freesteam_region2_h_pT(S.R2.p,S.R2.T);
        x=1;

        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        //Cl: note: in FreeStream, kappa=1/V*(dV/dP)_T=const is called kappaT (in this region)
        kappa=freesteam_region2_kappaT_pT(S.R2.p,S.R2.T);
        beta=freesteam_region2_alphav_pT(S.R2.p,S.R2.T);
        cp=freesteam_region2_cp_pT(S.R2.p,S.R2.T);

        //CL: getting derivatives using Bridgmans table
        //CL: psi=(drho/dp)_h=const
        //CL: drhodh=(drho/dh)_p=const
        psi=-((T*beta*beta-beta)/cp-kappa*rho);
        drhodh=-rho*beta/cp;

        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
        lambda=freesteam_k_rhoT(rho,T);
        alpha=lambda/cp; //Cl: Important info -->alpha= thermal diffusivity time density
    }
    //CL: supercritial fluid
    else if (region==3)
    {
        scalar gamma,cv;

        rho=S.R3.rho;
        T=S.R3.T;
        p=freesteam_region3_p_rhoT(S.R3.rho,S.R3.T);
        h=freesteam_region3_h_rhoT(S.R3.rho,S.R3.T);

        //CL= when h<h @ critical point -->x=0 else x=1
        if (h<2084256.263)
        {
            x=0;
        }
        else
        {
            x=1;
        }

        //Cl: note: beta=1/V*(dV/dP)_P=const
        //Cl: note: kappa=1/V*(dV/dP)_T=const
        //Cl: note: in FreeStream, gamma=1/p*(dp/dT)_v=const is called alphap (in this region)
        gamma=freesteam_region3_alphap_rhoT(S.R3.rho,S.R3.T);
        cp=freesteam_region3_cp_rhoT(S.R3.rho,S.R3.T);
        cv=freesteam_region3_cv_rhoT(S.R3.rho,S.R3.T);
        beta=(cp-cv)/(S.R3.T/S.R3.rho*p*gamma);
        kappa=(cp-cv)/(S.R3.T/S.R3.rho*p*p*gamma*gamma);

        //CL: getting derivatives using Bridgmans table
        //CL: psi=(drho/dp)_h=const
        //CL: drhodh=(drho/dh)_p=const
        psi=-((T*beta*beta-beta)/cp-kappa*rho);
        drhodh=-rho*beta/cp;


        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
        lambda=freesteam_k_rhoT(rho,T);
        alpha=lambda/cp; //Cl: Important info -->alpha= thermal diffusivity time density

    }
    //inside the vapor dome
    else if (region==4)
    {
        scalar rhov,rhol,betav,betal,kappav,kappal,vv,vl,cpl,cpv,hl,hv,cp;
        scalar dvldp,dvvdp,dhldp,dhvdp;
        scalar dpdT,dvdh,dvdp,dxdp;

        SteamState Sl,Sv;

        x=S.R4.x;
        T=S.R4.T;
        rho=1/freesteam_region4_v_Tx(S.R4.T,S.R4.x);
        h=freesteam_region4_h_Tx(S.R4.T,S.R4.x);
        p=freesteam_region4_psat_T(S.R4.T);
        cp=freesteam_region4_cp_Tx(S.R4.T,S.R4.x);


        //CL: Getting density on the vapour and liquid lines
        rhov=freesteam_region4_rhog_T(S.R4.T);
        rhol=freesteam_region4_rhof_T(S.R4.T);
        vv=1/rhov;
        vl=1/rhol;

        //CL: getting derivatives --> this is a bit tricky inside the vapor dome

        dpdT=freesteam_region4_dpsatdT_T(S.R4.T);

        // getting the states outside the vapour dome
        Sl=freesteam_set_pv(p,vl-0.0000001);  //inside region 1
        Sv=freesteam_set_pv(p,vv+0.0000001);  //inside region 2

        kappal=freesteam_region1_kappaT_pT(Sl.R1.p,Sl.R1.T);
        kappav=freesteam_region2_kappaT_pT(Sv.R2.p,Sv.R2.T);

        betal=freesteam_region1_alphav_pT(Sl.R1.p,Sl.R1.T);
        betav=freesteam_region2_alphav_pT(Sv.R2.p,Sv.R2.T);

        cpl=freesteam_region1_cp_pT(Sl.R1.p,Sl.R1.T);
        cpv=freesteam_region2_cp_pT(Sv.R2.p,Sv.R2.T);

        hl=freesteam_region1_h_pT(Sl.R1.p,Sl.R1.T);
        hv=freesteam_region2_h_pT(Sv.R2.p,Sv.R2.T);


        //calculation derviatives on liquid and vapour line
        dvldp=betal*vl/dpdT-kappal*vl;
        dvvdp=betav*vv/dpdT-kappav*vv;

        dhldp=vl*(1-betal*Sl.R1.T)+cpl/dpdT;
        dhvdp=vv*(1-betav*Sv.R2.T)+cpv/dpdT;

        dxdp=-dhldp/(hv-hl)
                 +(h-hl)/((hv-hl)*(hv-hl))
                     *(dhvdp-dhldp);

        //CL: psi=(drho/dp)_h=const
        dvdp=dvldp+(dvvdp-dvldp)*x+(vv-vl)*dxdp;
        psi=-rho*rho*dvdp;

        //CL: drhodh=(drho/dh)_p=const
        dvdh=(vv-vl)/(hv-hl);
        drhodh=-rho*rho*dvdh;

        //CL: getting transport properties
        mu=freesteam_mu_rhoT(rho, T);
        lambda=freesteam_k_rhoT(rho,T);
        alpha=lambda/cp; //Cl: Important info -->alpha= thermal diffusivity time density
    }
    else
    {
        std::cout<<"IAPWS-IF97 error, outside the regions 1-4"<<std::endl;
    }
}


//CL: returns density for given pressure and temperature
Foam::scalar Foam::rho_pT(scalar p,scalar T)
{
    return 1/freesteam_v(freesteam_set_pT(p,T));
}


//CL: returns density for given pressure and enthalpy
Foam::scalar Foam::rho_ph(scalar p,scalar h)
{
    return 1/freesteam_v(freesteam_set_ph(p,h));
}


//CL: returns Cp(heat capacity @ contant pressure) for given pressure and temperature
Foam::scalar Foam::cp_pT(scalar p,scalar T)
{
    return freesteam_cp(freesteam_set_pT(p,T));
}


//CL: returns Cp(heat capacity @ contant pressure) for given pressure and enthalpy
Foam::scalar Foam::cp_ph(scalar p,scalar h)
{
    return freesteam_cp(freesteam_set_ph(p,h));
}

//CL: returns Cv (heat capacity @ contant volume) for given pressure and temperature
Foam::scalar Foam::cv_pT(scalar p,scalar T)
{
    return freesteam_cv(freesteam_set_pT(p,T));
}


//CL: returns Cv (heat capacity @ contant volume) for given pressure and enthalpy
Foam::scalar Foam::cv_ph(scalar p,scalar h)
{
    return freesteam_cv(freesteam_set_ph(p,h));
}


//CL: returns enthalpy for given pressure and temperature
Foam::scalar Foam::h_pT(scalar p,scalar T)
{
    return freesteam_h(freesteam_set_pT(p,T));
}


//CL: returns temperature for given pressure and enthalpy
Foam::scalar Foam::T_ph(scalar p,scalar h)
{
    return freesteam_T(freesteam_set_ph(p,h));
}


//CL: psiH=(drho/dp)_h=const
Foam::scalar Foam::psiH_pT(scalar p,scalar T)
{
    return psiH(freesteam_set_pT(p,T));
}


//CL: psiH=(drho/dp)_h=const
Foam::scalar Foam::psiH_ph(scalar p,scalar h)
{
    return psiH(freesteam_set_ph(p,h));
}


//CL: psiH=(drho/dp)_h=const
Foam::scalar Foam::psiH(SteamState S)
{
    label region;
    scalar kappa,cp,beta,psiH,rho;

    region=freesteam_region(S);

    //CL:liquid phase
    if (region==1)
    {
        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        //Cl: note: in FreeStream, kappa=1/V*(dV/dP)_T=const is called kappaT (in this region)
        kappa=freesteam_region1_kappaT_pT(S.R1.p,S.R1.T);
        beta=freesteam_region1_alphav_pT(S.R1.p,S.R1.T);
        cp=freesteam_region1_cp_pT(S.R1.p,S.R1.T);
        rho=1/freesteam_region1_v_pT(S.R1.p,S.R1.T);

        //CL: getting derivatives using Bridgmans table
        //CL: psiH=(drho/dp)_h=const
        psiH=-((S.R1.T*beta*beta-beta)/cp-kappa*rho);
    }
    //CL:vapor phase
    else if (region==2)
    {
        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        //Cl: note: in FreeStream, kappa=1/V*(dV/dP)_T=const is called kappaT (in this region)
        kappa=freesteam_region2_kappaT_pT(S.R2.p,S.R2.T);
        beta=freesteam_region2_alphav_pT(S.R2.p,S.R2.T);
        cp=freesteam_region2_cp_pT(S.R2.p,S.R2.T);
        rho=1/freesteam_region2_v_pT(S.R2.p,S.R2.T);

        //CL: getting derivatives using Bridgmans table
        //CL: psiH=(drho/dp)_h=const
        psiH=-((S.R2.T*beta*beta-beta)/cp-kappa*rho);
    }
    //CL:supercritical fluid
    else if (region==3)
    {

        scalar gamma,cv,p;

        rho=S.R3.rho;
        p=freesteam_region3_p_rhoT(S.R3.rho,S.R3.T);

        //Cl: note: beta=1/V*(dV/dP)_P=const
        //Cl: note: kappa=1/V*(dV/dP)_T=const
        //Cl: note: in FreeStream, gamma=1/p*(dp/dT)_v=const is called alphap (in this region
        gamma=freesteam_region3_alphap_rhoT(S.R3.rho,S.R3.T);
        cp=freesteam_region3_cp_rhoT(S.R3.rho,S.R3.T);
        cv=freesteam_region3_cv_rhoT(S.R3.rho,S.R3.T);
        beta=(cp-cv)/(S.R3.T/S.R3.rho*p*gamma);
        kappa=(cp-cv)/(S.R3.T/S.R3.rho*p*p*gamma*gamma);

        //CL: getting derivatives using Bridgmans table
        //CL: psiH=(drho/dp)_h=const
        psiH=-((S.R3.T*beta*beta-beta)/cp-kappa*rho);

    }
    //inside the vapor dome
    else if (region==4)
    {
        scalar rhov,rhol,betav,betal,kappav,kappal,vv,vl,cpl,cpv,hl,hv,h,p;
        scalar dvldp,dvvdp,dhldp,dhvdp;
        scalar dpdT,dvdp,dxdp;

        SteamState Sl,Sv;

        rho=1/freesteam_region4_v_Tx(S.R4.T,S.R4.x);
        h=freesteam_region4_h_Tx(S.R4.T,S.R4.x);
        p=freesteam_region4_psat_T(S.R4.T);

        //CL: Getting density on the vapour and liquid lines
        rhov=freesteam_region4_rhog_T(S.R4.T);
        rhol=freesteam_region4_rhof_T(S.R4.T);
        vv=1/rhov;
        vl=1/rhol;

        //CL: getting derivatives --> this is a bit tricky in the vapor dome

        dpdT=freesteam_region4_dpsatdT_T(S.R4.T);

        // getting the states outside the vapour dome
        Sl=freesteam_set_pv(p,vl-0.0000001);  //inside region 1
        Sv=freesteam_set_pv(p,vv+0.0000001);  //inside region 2

        kappal=freesteam_region1_kappaT_pT(Sl.R1.p,Sl.R1.T);
        kappav=freesteam_region2_kappaT_pT(Sv.R2.p,Sv.R2.T);

        betal=freesteam_region1_alphav_pT(Sl.R1.p,Sl.R1.T);
        betav=freesteam_region2_alphav_pT(Sv.R2.p,Sv.R2.T);

        cpl=freesteam_region1_cp_pT(Sl.R1.p,Sl.R1.T);
        cpv=freesteam_region2_cp_pT(Sv.R2.p,Sv.R2.T);

        hl=freesteam_region1_h_pT(Sl.R1.p,Sl.R1.T);
        hv=freesteam_region2_h_pT(Sv.R2.p,Sv.R2.T);

        //calculation derviatives on liquid and vapour line
        dvldp=betal*vl/dpdT-kappal*vl;
        dvvdp=betav*vv/dpdT-kappav*vv;

        dhldp=vl*(1-betal*Sl.R1.T)+cpl/dpdT;
        dhvdp=vv*(1-betav*Sv.R2.T)+cpv/dpdT;

        dxdp=-dhldp/(hv-hl)
                 +(h-hl)/((hv-hl)*(hv-hl))
                     *(dhvdp-dhldp);

        //CL: psiH=(drho/dp)_h=const
        dvdp=dvldp+(dvvdp-dvldp)*S.R4.x+(vv-vl)*dxdp;
        psiH=-rho*rho*dvdp;
    }
    else
    {
        Info<<"IAPWS-IF97.C error, outside the regions 1-4"<<endl;

        //Keep compiler happy
        psiH = 0;
    }

    return psiH;
}


//CL: drhodh=(drho/dh)_p=const
Foam::scalar Foam::drhodh_pT(scalar p,scalar T)
{
    return drhodh(freesteam_set_pT(p,T));
}


//CL: drhodh=(drho/dh)_p=const
Foam::scalar Foam::drhodh_ph(scalar p,scalar h)
{
    return drhodh(freesteam_set_ph(p,h));
}


//CL: drhodh=(drho/dh)_p=const
Foam::scalar Foam::drhodh(SteamState S)
{
    label region;
    scalar cp,beta,drhodh,rho;

    region=freesteam_region(S);

    if (region==1)
    {
        rho=1/freesteam_region1_v_pT(S.R1.p,S.R1.T);

        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        beta=freesteam_region1_alphav_pT(S.R1.p,S.R1.T);
        cp=freesteam_region1_cp_pT(S.R1.p,S.R1.T);

        //CL: getting derivatives using Bridgmans table
        //CL: drhodh=(drho/dh)_p=const
        drhodh=-rho*beta/cp;
    }
    else if (region==2)
    {
        rho=1/freesteam_region2_v_pT(S.R2.p,S.R2.T);

        //Cl: note: in FreeStream, beta=1/V*(dV/dP)_P=const is called alphaV (in this region)
        //Cl: note: in FreeStream, kappa=1/V*(dV/dP)_T=const is called kappaT (in this region)
        beta=freesteam_region2_alphav_pT(S.R2.p,S.R2.T);
        cp=freesteam_region2_cp_pT(S.R2.p,S.R2.T);

        //CL: getting derivatives using Bridgmans table
        //CL: drhodh=(drho/dh)_p=const
        drhodh=-rho*beta/cp;
    }
    else if (region==3)
    {
        scalar gamma,cv,p;

        p=freesteam_region3_p_rhoT(S.R3.rho,S.R3.T);

        //Cl: note: beta=1/V*(dV/dP)_P=const
        //Cl: note: in FreeStream, gamma=1/p*(dp/dT)_v=const is called alphap (in this region
        gamma=freesteam_region3_alphap_rhoT(S.R3.rho,S.R3.T);
        cp=freesteam_region3_cp_rhoT(S.R3.rho,S.R3.T);
        cv=freesteam_region3_cv_rhoT(S.R3.rho,S.R3.T);
        beta=(cp-cv)/(S.R3.T/S.R3.rho*p*gamma);

        //CL: getting derivatives using Bridgmans table
        //CL: drhodh=(drho/dh)_p=const
        drhodh=-S.R3.rho*beta/cp;
    }
    else if (region==4)
    {
        scalar vv,vl,hl,hv,p;
        SteamState Sl,Sv;

        rho=1/freesteam_region4_v_Tx(S.R4.T,S.R4.x);
        p=freesteam_region4_psat_T(S.R4.T);

        //CL: Getting density on the vapour and liquid lines
        vv=1/freesteam_region4_rhog_T(S.R4.T);
        vl=1/freesteam_region4_rhof_T(S.R4.T);

        // getting the states outside the vapour dome
        Sl=freesteam_set_pv(p,vl-0.0000001);  //inside region 1
        Sv=freesteam_set_pv(p,vv+0.0000001);  //inside region 2

        hl=freesteam_region1_h_pT(Sl.R1.p,Sl.R1.T);
        hv=freesteam_region2_h_pT(Sv.R2.p,Sv.R2.T);

        //CL: drhodh=(drho/dh)_p=const
        drhodh=-rho*rho*(vv-vl)/(hv-hl);
    }
    else
    {
         Info<<"IAPWS-IF97.C error, outside the regions 1-4"<<endl;

         // Keep compiler happy
         drhodh = 0;
    }

    return drhodh;
}


