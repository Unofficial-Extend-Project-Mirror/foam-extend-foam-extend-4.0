// Parametrized test case for the ERCOFTAC diffuser.

// Created by Omar Bounous 

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926536)
m4_define(sqr,0.707106781186548)

//Geometry
m4_define(openingAngle, 10.0)
m4_define(diffuserLength, 0.51)
m4_define(extensionLength, 0.59)
m4_define(rIn,0.13)
m4_define(rOut, calc(rIn+diffuserLength*tan(openingAngle*pi/180.0)))

//Grid points (integers!):
m4_define(rNumberOfCells, 25)
m4_define(tNumberOfCells, 20)
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
//Plane A:
(0 0 zA) vlabel(A0)
(rA 0 zA) vlabel(A1)
(0 rA zA) vlabel(A2)
(-rA 0 zA) vlabel(A3)
(0 -rA zA) vlabel(A4)

//Plane B:
(0 0 zB) vlabel(B0)
(rB 0 zB) vlabel(B1)
(0 rB zB) vlabel(B2)
(-rB 0 zB) vlabel(B3)
(0 -rB zB) vlabel(B4)

//Plane C:
(0 0 zC) vlabel(C0)
(rC 0 zC) vlabel(C1)
(0 rC zC) vlabel(C2)
(-rC 0 zC) vlabel(C3)
(0 -rC zC) vlabel(C4)

//Plane D:
(0 0 zD) vlabel(D0)
(rD 0 zD) vlabel(D1)
(0 rD zD) vlabel(D2)
(-rD 0 zD) vlabel(D3)
(0 -rD zD) vlabel(D4)

//Plane E:
(0 0 zE) vlabel(E0)
(rE 0 zE) vlabel(E1)
(0 rE zE) vlabel(E2)
(-rE 0 zE) vlabel(E3)
(0 -rE zE) vlabel(E4)

//Plane F:
(0 0 zF) vlabel(F0)
(rF 0 zF) vlabel(F1)
(0 rF zF) vlabel(F2)
(-rF 0 zF) vlabel(F3)
(0 -rF zF) vlabel(F4)

);

// Defining blocks:
blocks
(
    //Blocks between plane A and plane B:
    // block0 - positive x and y 
    hex (A0 A1 A2 A0 B0 B1 B2 B0) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block1 - negative x positive y 
    hex (A0 A2 A3 A0 B0 B2 B3 B0) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block2 - negative x and  y 
    hex (A0 A3 A4 A0 B0 B3 B4 B0) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block3 - positive x negative y 
    hex (A0 A4 A1 A0 B0 B4 B1 B0) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane B and plane C:
    // block0 - positive x and y 
    hex (B0 B1 B2 B0 C0 C1 C2 C0) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block1 - negative x positive y 
    hex (B0 B2 B3 B0 C0 C2 C3 C0) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block2 - negative x and  y 
    hex (B0 B3 B4 B0 C0 C3 C4 C0) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block3 - positive x negative y 
    hex (B0 B4 B1 B0 C0 C4 C1 C0) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane C and plane D:
    // block0 - positive x and y 
    hex (C0 C1 C2 C0 D0 D1 D2 D0) CD
    (rNumberOfCells tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block1 - negative x positive y 
    hex (C0 C2 C3 C0 D0 D2 D3 D0) CD
    (rNumberOfCells tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block2 - negative x and  y 
    hex (C0 C3 C4 C0 D0 D3 D4 D0) CD
    (rNumberOfCells tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block3 - positive x negative y 
    hex (C0 C4 C1 C0 D0 D4 D1 D0) CD
    (rNumberOfCells tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane D and plane E:
    // block0 - positive x and y 
    hex (D0 D1 D2 D0 E0 E1 E2 E0) DE
    (rNumberOfCells tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block1 - negative x positive y 
    hex (D0 D2 D3 D0 E0 E2 E3 E0) DE
    (rNumberOfCells tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block2 - negative x and  y 
    hex (D0 D3 D4 D0 E0 E3 E4 E0) DE
    (rNumberOfCells tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block3 - positive x negative y 
    hex (D0 D4 D1 D0 E0 E4 E1 E0) DE
    (rNumberOfCells tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading 1 1)

    //Blocks between plane E and plane F:
    // block0 - positive x and y 
    hex (E0 E1 E2 E0 F0 F1 F2 F0) EF
    (rNumberOfCells tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block1 - negative x positive y 
    hex (E0 E2 E3 E0 F0 F2 F3 F0) EF
    (rNumberOfCells tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block2 - negative x and  y 
    hex (E0 E3 E4 E0 F0 F3 F4 F0) EF
    (rNumberOfCells tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block3 - positive x negative y 
    hex (E0 E4 E1 E0 F0 F4 F1 F0) EF
    (rNumberOfCells tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading 1 1)
);

edges
(
    //Plane A:
    arc A1 A2 (calc(rA*sqr) calc(rA*sqr) zA)
    arc A2 A3 (-calc(rA*sqr) calc(rA*sqr) zA)
    arc A3 A4 (-calc(rA*sqr) -calc(rA*sqr) zA)
    arc A4 A1 (calc(rA*sqr) -calc(rA*sqr) zA)

    //Plane B:
    arc B1 B2 (calc(rB*sqr) calc(rB*sqr) zB)
    arc B2 B3 (-calc(rB*sqr) calc(rB*sqr) zB)
    arc B3 B4 (-calc(rB*sqr) -calc(rB*sqr) zB)
    arc B4 B1 (calc(rB*sqr) -calc(rB*sqr) zB)

    //Plane C:
    arc C1 C2 (calc(rC*sqr) calc(rC*sqr) zC)
    arc C2 C3 (-calc(rC*sqr) calc(rC*sqr) zC)
    arc C3 C4 (-calc(rC*sqr) -calc(rC*sqr) zC)
    arc C4 C1 (calc(rC*sqr) -calc(rC*sqr) zC)

    //Plane D:
    arc D1 D2 (calc(rD*sqr) calc(rD*sqr) zD)
    arc D2 D3 (-calc(rD*sqr) calc(rD*sqr) zD)
    arc D3 D4 (-calc(rD*sqr) -calc(rD*sqr) zD)
    arc D4 D1 (calc(rD*sqr) -calc(rD*sqr) zD)

    //Plane E:
    arc E1 E2 (calc(rE*sqr) calc(rE*sqr) zE)
    arc E2 E3 (-calc(rE*sqr) calc(rE*sqr) zE)
    arc E3 E4 (-calc(rE*sqr) -calc(rE*sqr) zE)
    arc E4 E1 (calc(rE*sqr) -calc(rE*sqr) zE)

    //Plane F:
    arc F1 F2 (calc(rF*sqr) calc(rF*sqr) zF)
    arc F2 F3 (-calc(rF*sqr) calc(rF*sqr) zF)
    arc F3 F4 (-calc(rF*sqr) -calc(rF*sqr) zF)
    arc F4 F1 (calc(rF*sqr) -calc(rF*sqr) zF)

);

// Defining patches:
patches
(
    patch inlet
    (
       (A0 A2 A1 A0)
       (A0 A1 A4 A0)
       (A0 A4 A3 A0)
       (A0 A3 A2 A0)
    )
    patch outlet
    (
       (F0 F1 F2 F0)
       (F0 F2 F3 F0)
       (F0 F3 F4 F0)
       (F0 F4 F1 F0)
    )
    wall wallProlongation
    (
      (E1 E2 F2 F1)
      (E2 E3 F3 F2)
      (E3 E4 F4 F3)
      (E4 E1 F1 F4)
    )
    wall wallDiffuser
    (
      (D1 D2 E2 E1)
      (D2 D3 E3 E2)
      (D3 D4 E4 E3)
      (D4 D1 E1 E4)
    )
    wall statSwirlWall
    (
      (B1 B2 C2 C1)
      (B2 B3 C3 C2)
      (B3 B4 C4 C3)
      (B4 B1 C1 C4)
      (C1 C2 D2 D1)
      (C2 C3 D3 D2)
      (C3 C4 D4 D3)
      (C4 C1 D1 D4)
    )
   wall rotSwirlWall
   (
      (A1 A2 B2 B1)
      (A2 A3 B3 B2)
      (A3 A4 B4 B3)
      (A4 A1 B1 B4)
    )
);

mergePatchPairs 
(
);

// ************************************************************************* //
