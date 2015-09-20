/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     3.2                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
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
// General macros to create 2D/extruded-2D meshes

changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
//define(calc, [esyscmd(echo $1 | bc | tr -d \\n)])
define(VCOUNT, 0)
define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])
define(pi, calc(3.14159265/20))

define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))
define(quad2D, ($1b $2b $2t $1t))
define(frontQuad, ($1t $2t $3t $4t))
define(backQuad, ($1b $4b $3b $2b))

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

// HUB AND SHROUD RADIUS
// Hub radius (m)
define(hr, 0.05)
// Shroud radius (m)
define(sr, 0.1)

// RUNNER REGION GEOMETRY AND MESH PROPERTIES
// Runner inlet axial length (m)
define(RUial, 0.02)
// Runner axial length (m)
define(RUal, 0.1)
// Runner outlet axial length (m)
define(RUoal, 0.02)
// Number of runner blades per 360 degrees (integer!)
define(RUnb, 5)
// Number of cells in radial direction in runner
define(RUrc, 10)
// Number of cells in tangential direction between runner blades
define(RUtc, 10)
// Number of cells in axial direction at runner inlet
define(RUiac, 2)
// Number of cells in axial direction between runner blades
define(RUbac, 10)
// Number of cells in axial direction at runner outlet
define(RUoac, 2)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TANGENTIAL PITCHES (RADIANS)
// Runner region
define(RUp, calc(2*pi/RUnb))

// TANGENTIAL SHIFTS BETWEEN AXIAL LEVELS (BOTTOM-UP)
// Runner region
// Tangential shift from level RU0 to RU1
define(RUts01, calc(-1/10*RUp))
// Tangential shift from level RU1 to RU2
define(RUts12, calc(-4/5*RUp))
// Tangential shift from level RU2 to RU3
define(RUts23, calc(-1/10*RUp))

// AXIAL/TANGENTIAL BASE POINTS FOR EACH LEVEL (BOTTOM-UP):
// (CENTER OF RUNNER SET TO THETA=0, Z=0)
// Runner:
define(RUa0, calc(-RUoal-0.5*RUal)) //Center runner
define(RUt0, calc(-0.5*RUp-(0.5*RUts12))) //Center runner
define(RUt1, calc(RUt0+RUts01))
define(RUt2, calc(RUt1+RUts12))
define(RUt3, calc(RUt2+RUts23))

vertices //(radial [m], tangential [radians], axial [m])
(
//Runner hub:
    (hr RUt0 RUa0) vlabel(RU0lb)
    (hr calc(RUt0+RUp) RUa0) vlabel(RU0rb)
    (hr RUt1 calc(RUa0+RUoal)) vlabel(RU1lb)
    (hr calc(RUt1+RUp) calc(RUa0+RUoal)) vlabel(RU1rb)
    (hr RUt2 calc(RUa0+RUoal+RUal)) vlabel(RU2lb)
    (hr calc(RUt2+RUp) calc(RUa0+RUoal+RUal)) vlabel(RU2rb)
    (hr RUt3 calc(RUa0+RUoal+RUal+RUial)) vlabel(RU3lb)
    (hr calc(RUt3+RUp) calc(RUa0+RUoal+RUal+RUial)) vlabel(RU3rb)

//Runner shroud:
    (sr RUt0 RUa0) vlabel(RU0lt)
    (sr calc(RUt0+RUp) RUa0) vlabel(RU0rt)
    (sr RUt1 calc(RUa0+RUoal)) vlabel(RU1lt)
    (sr calc(RUt1+RUp) calc(RUa0+RUoal)) vlabel(RU1rt)
    (sr RUt2 calc(RUa0+RUoal+RUal)) vlabel(RU2lt)
    (sr calc(RUt2+RUp) calc(RUa0+RUoal+RUal)) vlabel(RU2rt)
    (sr RUt3 calc(RUa0+RUoal+RUal+RUial)) vlabel(RU3lt)
    (sr calc(RUt3+RUp) calc(RUa0+RUoal+RUal+RUial)) vlabel(RU3rt)
);

blocks
(
//Runner:
    hex2D(RU0l, RU0r, RU1r, RU1l)
    rotor
    (RUtc RUoac RUrc)
    simpleGrading (1 1 1)

    hex2D(RU1l, RU1r, RU2r, RU2l)
    rotor
    (RUtc RUbac RUrc)
    simpleGrading (1 1 1)

    hex2D(RU2l, RU2r, RU3r, RU3l)
    rotor
    (RUtc RUiac RUrc)
    simpleGrading (1 1 1)
);

edges // Inappropriate with arc due to coordinate conversion
(
//Runner
    spline RU1lt RU2lt
    (
        (sr calc(RUt1+0.65*(RUt2-(RUt1))) calc(RUa0+RUoal+0.5*RUal))
    )
    spline RU1lb RU2lb
    (
        (hr calc(RUt1+0.65*(RUt2-(RUt1))) calc(RUa0+RUoal+0.5*RUal))
    )
    spline RU1rt RU2rt
    (
        (sr calc(RUt1+RUp+0.75*(RUt2-(RUt1))) calc(RUa0+RUoal+0.5*RUal))
    )
    spline RU1rb RU2rb
    (
        (hr calc(RUt1+RUp+0.75*(RUt2-(RUt1))) calc(RUa0+RUoal+0.5*RUal))
    )
);

boundary
(
    RUINLET
    {
        type            patch;
        faces
        (
            quad2D(RU3r, RU3l)
        );
    }

    RUOUTLET
    {
        type            patch;
        faces
        (
            quad2D(RU0l, RU0r)
        );
    }

    RUCYCLIC1
    {
        type             cyclicGgi;
        shadowPatch      RUCYCLIC2;
        zone             RUCYCLIC1Zone;
        bridgeOverlap    false;
        rotationAxis     (0 0 1);
        rotationAngle    72;
        separationOffset (0 0 0);
        faces
        (
            quad2D(RU1l, RU0l)
            quad2D(RU3l, RU2l)
        );
    }

    RUCYCLIC2
    {
        type             cyclicGgi;
        shadowPatch      RUCYCLIC1;
        zone             RUCYCLIC2Zone;
        bridgeOverlap    false;
        rotationAxis     (0 0 1);
        rotationAngle    -72;
        separationOffset (0 0 0);
        faces
        (
            quad2D(RU0r, RU1r)
            quad2D(RU2r, RU3r)
        );
    }

    RUBLADE
    {
        type wall;
        faces
        (
            quad2D(RU2l, RU1l)
            quad2D(RU1r, RU2r)
        );
    }

    RUHUB
    {
        type wall;
        faces
        (
            backQuad(RU0l, RU0r, RU1r, RU1l)
            backQuad(RU1l, RU1r, RU2r, RU2l)
            backQuad(RU2l, RU2r, RU3r, RU3l)
        );
    }

    RUSHROUD
    {
        type wall;
        faces
        (
            frontQuad(RU0l, RU0r, RU1r, RU1l)
            frontQuad(RU1l, RU1r, RU2r, RU2l)
            frontQuad(RU2l, RU2r, RU3r, RU3l)
        );
    }
);

// ************************************************************************* //
