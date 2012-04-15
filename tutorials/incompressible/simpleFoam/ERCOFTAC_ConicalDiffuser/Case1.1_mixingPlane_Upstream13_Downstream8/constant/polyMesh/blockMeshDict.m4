// Parametrized test case for the ERCOFTAC diffuser.

// Created by Omar Bounous 

// Modified by Martin Beaudoin (11/2009)

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])

m4_define(pi, 3.14159265358979323844)
m4_define(rad, [calc($1*pi/180.0)])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(nbrPassageUpstream, 13)
m4_define(nbrPassageDownstream, 8)

// Make sure we are properly aligned for sampling inside the mesh volume
//m4_define(startAngleOffsetA_BB, rad(calc(90.0 - 360.0/nbrPassageUpstream/2.0)))
m4_define(startAngleOffsetB_F,  rad(calc(90.0 - 360.0/nbrPassageDownstream/2.0)))
m4_define(startAngleOffsetA_BB, rad(0.0))
//m4_define(startAngleOffsetB_F,  rad(5.0))


// Angle span for section A_BB
m4_define(startAngleSectionA_BB, calc(startAngleOffsetA_BB + rad(0.0)))
m4_define(stopAngleSectionA_BB,  calc(startAngleOffsetA_BB + rad(calc(360/nbrPassageUpstream))))
m4_define(angleSpanSectionA_BB,  calc(stopAngleSectionA_BB - startAngleSectionA_BB))

// Angle span for section B_F
m4_define(startAngleSectionB_F, calc(startAngleOffsetB_F + rad(0.0)))
m4_define(stopAngleSectionB_F,  calc(startAngleOffsetB_F + rad(calc(360/nbrPassageDownstream))))
m4_define(angleSpanSectionB_F,  calc(stopAngleSectionB_F - startAngleSectionB_F))


//Geometry
m4_define(openingAngle, 10.0)
m4_define(diffuserLength, 0.51)
m4_define(extensionLength, 0.59)
m4_define(rIn,0.13)
m4_define(rOut, calc(rIn+diffuserLength*tan(openingAngle*pi/180.0)))

//Grid points (integers!):
// For a better mesh resolution in the radial and tangential direction, 
// play with rNumberOfCells and tNumberOfCells
m4_define(rUpstreamNumberOfCells, 21)
m4_define(rDownstreamNumberOfCells, 27)
m4_define(tNumberOfCells, 20)
//m4_define(rNumberOfCells, 3)
//m4_define(tNumberOfCells, 4)
//
m4_define(zABnumberOfCells, 10)
m4_define(zBCnumberOfCells, 4)
m4_define(zCDnumberOfCells, 6)
m4_define(zDEnumberOfCells, 30)
m4_define(zEFnumberOfCells, 10)
m4_define(rGrading, 0.2)

//Plane A:
m4_define(zA, -0.50)
m4_define(rA, rIn)

//Plane B:
m4_define(zB, -0.10)
m4_define(rB, rIn)

//Plane BB:     // This is where we put the mixingPlane interface
m4_define(zBB, -0.10)
m4_define(rBB, rIn)

//Plane C:
m4_define(zC, -0.025)
m4_define(rC, rIn)

//Plane D:
m4_define(zD, 0)
m4_define(rD, rIn)

//Plane E:
m4_define(zE, diffuserLength)
m4_define(rE, rOut)

//Plane F:
m4_define(zF, calc(diffuserLength+extensionLength))
m4_define(rF, rOut)

/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.5                                   |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
//Plane A: Only a sector
(                                0                                   0   zA)  vlabel(A0)
(calc(rA*cos(startAngleSectionA_BB)) calc(rA*sin(startAngleSectionA_BB)) zA)  vlabel(A1)
(calc(rA*cos( stopAngleSectionA_BB)) calc(rA*sin( stopAngleSectionA_BB)) zA)  vlabel(A2)

//Plane BB: Only a sector
(                                 0                                    0   zBB)  vlabel(BB0)
(calc(rBB*cos(startAngleSectionA_BB)) calc(rBB*sin(startAngleSectionA_BB)) zBB)  vlabel(BB1)
(calc(rBB*cos( stopAngleSectionA_BB)) calc(rBB*sin( stopAngleSectionA_BB)) zBB)  vlabel(BB2)

//Plane B: Only a sector
//(  0   0 zB) vlabel(B0)
//( rB   0 zB) vlabel(B1)
//(  0  rB zB) vlabel(B2)
//(-rB   0 zB) vlabel(B3)
//(  0 -rB zB) vlabel(B4)
(                               0                                  0   zB)  vlabel(B0)
(calc(rB*cos(startAngleSectionB_F)) calc(rB*sin(startAngleSectionB_F)) zB)  vlabel(B1)
(calc(rB*cos( stopAngleSectionB_F)) calc(rB*sin( stopAngleSectionB_F)) zB)  vlabel(B2)

//Plane C: Only a sector
//(0 0 zC) vlabel(C0)
//(rC 0 zC) vlabel(C1)
//(0 rC zC) vlabel(C2)
//(-rC 0 zC) vlabel(C3)
//(0 -rC zC) vlabel(C4)
(                               0                                  0   zC)  vlabel(C0)
(calc(rC*cos(startAngleSectionB_F)) calc(rC*sin(startAngleSectionB_F)) zC)  vlabel(C1)
(calc(rC*cos( stopAngleSectionB_F)) calc(rC*sin( stopAngleSectionB_F)) zC)  vlabel(C2)

//Plane D:
//(0 0 zD) vlabel(D0)
//(rD 0 zD) vlabel(D1)
//(0 rD zD) vlabel(D2)
//(-rD 0 zD) vlabel(D3)
//(0 -rD zD) vlabel(D4)
(                               0                                  0   zD)  vlabel(D0)
(calc(rD*cos(startAngleSectionB_F)) calc(rD*sin(startAngleSectionB_F)) zD)  vlabel(D1)
(calc(rD*cos( stopAngleSectionB_F)) calc(rD*sin( stopAngleSectionB_F)) zD)  vlabel(D2)

//Plane E:
//(0 0 zE) vlabel(E0)
//(rE 0 zE) vlabel(E1)
//(0 rE zE) vlabel(E2)
//(-rE 0 zE) vlabel(E3)
//(0 -rE zE) vlabel(E4)
(                               0                                  0   zE)  vlabel(E0)
(calc(rE*cos(startAngleSectionB_F)) calc(rE*sin(startAngleSectionB_F)) zE)  vlabel(E1)
(calc(rE*cos( stopAngleSectionB_F)) calc(rE*sin( stopAngleSectionB_F)) zE)  vlabel(E2)

//Plane F:
//(0 0 zF) vlabel(F0)
//(rF 0 zF) vlabel(F1)
//(0 rF zF) vlabel(F2)
//(-rF 0 zF) vlabel(F3)
//(0 -rF zF) vlabel(F4)
(                               0                                  0   zF)  vlabel(F0)
(calc(rF*cos(startAngleSectionB_F)) calc(rF*sin(startAngleSectionB_F)) zF)  vlabel(F1)
(calc(rF*cos( stopAngleSectionB_F)) calc(rF*sin( stopAngleSectionB_F)) zF)  vlabel(F2)

);

// Defining blocks:
blocks
(
    //Blocks between plane A and plane BB:
    // block0 - positive x and y 
    hex (A0 A1 A2 A0 BB0 BB1 BB2 BB0) ABB
    (rUpstreamNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane B and plane C:
    // block0 - positive x and y 
    hex (B0 B1 B2 B0 C0 C1 C2 C0) BC
    (rDownstreamNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane C and plane D:
    // block0 - positive x and y 
    hex (C0 C1 C2 C0 D0 D1 D2 D0) CD
    (rDownstreamNumberOfCells tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane D and plane E:
    // block0 - positive x and y 
    hex (D0 D1 D2 D0 E0 E1 E2 E0) DE
    (rDownstreamNumberOfCells tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane E and plane F:
    // block0 - positive x and y 
    hex (E0 E1 E2 E0 F0 F1 F2 F0) EF
    (rDownstreamNumberOfCells tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading 1 1)
);

//vstartAngleOffsetA_BB = startAngleOffsetA_BB
//vstartAngleSectionA_BB = startAngleSectionA_BB
//vstopAngleSectionA_BB = stopAngleSectionA_BB
//vangleSpanSectionA_BB = angleSpanSectionA_BB
//
//vstartAngleOffsetB_F = startAngleOffsetB_F
//vstartAngleSectionB_F = startAngleSectionB_F
//vstopAngleSectionB_F = stopAngleSectionB_F
//vangleSpanSectionB_F = angleSpanSectionB_F
//vCenterA_BB = A1 A2 (calc(rA*cos(startAngleOffsetA_BB + (angleSpanSectionA_BB)/2)) calc(rA*sin(startAngleOffsetA_BB + (angleSpanSectionA_BB)/2)) zA) 

edges 
(
    //Plane A:
    arc A1 A2 (calc(rA*cos(startAngleOffsetA_BB + (angleSpanSectionA_BB)/2)) calc(rA*sin(startAngleOffsetA_BB + (angleSpanSectionA_BB)/2)) zA) 

    //Plane BB:
    arc BB1 BB2 (calc(rBB*cos(startAngleOffsetA_BB + (angleSpanSectionA_BB)/2)) calc(rBB*sin(startAngleOffsetA_BB + (angleSpanSectionA_BB)/2)) zBB)

    //Plane B:
    arc B1 B2 (calc(rB*cos(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) calc(rB*sin(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) zB)

    //Plane C:
    arc C1 C2 (calc(rC*cos(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) calc(rC*sin(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) zC)

    //Plane D:
    arc D1 D2 (calc(rD*cos(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) calc(rD*sin(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) zD)

    //Plane E:
    arc E1 E2 (calc(rE*cos(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) calc(rE*sin(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) zE)

    //Plane F:
    arc F1 F2 (calc(rF*cos(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) calc(rF*sin(startAngleOffsetB_F + (angleSpanSectionB_F)/2)) zF)
);

// Defining patches:
patches
(
    patch inlet
    (
       (A0 A2 A1 A0)
    )
    patch outlet
    (
       (F0 F1 F2 F0)
    )
    wall wallProlongation
    (
      (E1 E2 F2 F1)
    )
    wall wallDiffuser
    (
      (D1 D2 E2 E1)
    )
    wall statSwirlWallB_C
    (
      (B1 B2 C2 C1)
      (C1 C2 D2 D1)
    )
   wall rotSwirlWallA_BB
   (
      (A1 A2 BB2 BB1)
   )

    cyclicGgi sideWallA_BB_cyclic1
    (
      (A0 BB0 BB1 A1)
    )

    cyclicGgi sideWallA_BB_cyclic2
    (
      (BB0 A0 A2 BB2)
    )

    mixingPlane B_UPSTREAM   // BB : master
    (
       (BB0 BB2 BB1 BB0)
    )

    mixingPlane B_DOWNSTREAM  // B : master
    (
       (B0 B2 B1 B0)
    )

    cyclicGgi sideWallB_C_cyclic1
    (
      (B0 C0 C1 B1)
    )

    cyclicGgi sideWallB_C_cyclic2
    (
      (C0 B0 B2 C2)
    )

    cyclicGgi sideWallC_D_cyclic1
    (
      (C0 D0 D1 C1)
    )

    cyclicGgi sideWallC_D_cyclic2
    (
      (D0 C0 C2 D2)
    )

    cyclicGgi sideWallD_E_cyclic1
    (
      (D0 E0 E1 D1)
    )

    cyclicGgi sideWallD_E_cyclic2
    (
      (E0 D0 D2 E2)
    )

    cyclicGgi sideWallE_F_cyclic1
    (
      (E0 F0 F1 E1)
    )

    cyclicGgi sideWallE_F_cyclic2
    (
      (F0 E0 E2 F2)
    )

);

mergePatchPairs 
(
);

// ************************************************************************* //
