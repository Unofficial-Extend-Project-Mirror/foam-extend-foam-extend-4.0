/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:  3.0                                |
|   \\  /    A nd           | Web:         http://www.extend-project.de       |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
//Plane A:
//Bottom of curved block
(calc(r1*cos(angleB)) calc(r1*sin(angleB)) zA) vlabel(A0)
(calc(r2*cos(angleB)) calc(r2*sin(angleB)) zA) vlabel(A1)
(calc(r1*cos(angleD)) calc(r1*sin(angleD)) zA) vlabel(A2)
(calc(r2*cos(angleD)) calc(r2*sin(angleD)) zA) vlabel(A3)

//Plane B:
//Top of curved block
(calc(r1*cos(angleB)) calc(r1*sin(angleB)) zB) vlabel(B0)
(calc(r2*cos(angleB)) calc(r2*sin(angleB)) zB) vlabel(B1)
(calc(r1*cos(angleD)) calc(r1*sin(angleD)) zB) vlabel(B2)
(calc(r2*cos(angleD)) calc(r2*sin(angleD)) zB) vlabel(B3)

//Plane A: Bottom of straight block
(calc(r2*cos(angleA)) calc(r2*sin(angleA)) zA) vlabel(B4)
(calc(r3*cos(angleA)) calc(r3*sin(angleA)) zA) vlabel(B5)
(calc(r2*cos(angleC)) calc(r2*sin(angleC)) zA) vlabel(B6)
(calc(r3*cos(angleC)) calc(r3*sin(angleC)) zA) vlabel(B7)

//Plane B: Top of straight block
(calc(r2*cos(angleA)) calc(r2*sin(angleA)) zB) vlabel(C0)
(calc(r3*cos(angleA)) calc(r3*sin(angleA)) zB) vlabel(C1)
(calc(r2*cos(angleC)) calc(r2*sin(angleC)) zB) vlabel(C2)
(calc(r3*cos(angleC)) calc(r3*sin(angleC)) zB) vlabel(C3)

);

blocks
(
    hex ( A0 A1 A3 A2 B0 B1 B3 B2 ) (BLOCKSIZE_UPSTREAM)   simpleGrading (1 1 grading)
    hex ( B4 B5 B7 B6 C0 C1 C3 C2 ) (BLOCKSIZE_DOWNSTREAM) simpleGrading (1 1 grading)
);

edges
(
    // --- PLANE A: Bottom of curved block
    arc  A0 A2  (calc(r1*cos((angleB+angleD)/2)) calc(r1*sin((angleB+angleD)/2)) zA)
    arc  A1 A3  (calc(r2*cos((angleB+angleD)/2)) calc(r2*sin((angleB+angleD)/2)) zA)

    // --- PLANE B: Top of curved block
    arc  B0 B2  (calc(r1*cos((angleB+angleD)/2)) calc(r1*sin((angleB+angleD)/2))  zB)
    arc  B1 B3  (calc(r2*cos((angleB+angleD)/2)) calc(r2*sin((angleB+angleD)/2)) zB)

    // --- PLANE A: Bottom of straight block
    arc  B4 B6  (calc(r2*cos((angleA+angleC)/2)) calc(r2*sin((angleA+angleC)/2)) zA)
    arc  B5 B7  (calc(r3*cos((angleA+angleC)/2)) calc(r3*sin((angleA+angleC)/2)) zA)

    // --- PLANE B: Top of straight block
    arc  C0 C2  (calc(r2*cos((angleA+angleC)/2)) calc(r2*sin((angleA+angleC)/2)) zB)
    arc  C1 C3  (calc(r3*cos((angleA+angleC)/2)) calc(r3*sin((angleA+angleC)/2)) zB)
);

patches
(
    patch outflow
    (
        ( B5 B7 C3 C1 )
    )
    patch inflow
    (
        ( A0 B0 B2 A2 )
    )

    cyclicGgi upstreamPerio1
    (
        (B4 B5 C1 C0)
    )

    cyclicGgi upstreamPerio2
    (
        (B6 C2 C3 B7)
    )

    mixingPlane upstreamMixingPlanePatch
    (
        ( B4 C0 C2 B6 )
    )

    mixingPlane downstreamMixingPlanePatch
    (
        ( A1 A3 B3 B1 )
    )

    symmetryPlane downstreamWall
    (
        ( A0 A2 A3 A1 )
        ( B0 B1 B3 B2 )
    )

    symmetryPlane upstreamWall
    (
        ( C1 C3 C2 C0)
        (B4 B6 B7 B5)
    )

    cyclicGgi downstreamPerio1
    (
        (A0 A1 B1 B0)
    )

    cyclicGgi downstreamPerio2
    (
        (A2 B2 B3 A3)
    )

);

mergePatchPairs
(
);
