/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     4.0                                |
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

// GUIDE VANE REGION GEOMETRY AND MESH PROPERTIES
// Guide vane inlet axial length (m)
define(GVial, 0.1)
// Guide vane axial length (m)
define(GVbal, 0.1)
// Guide vane outlet axial length (m)
define(GVoal, 0.02)
// Number of guide vanes per 360 degrees (integer!)
define(GVnb, 5)
// Number of cells in radial direction at guide vane
define(GVrc, 10)
// Number of cells in tangential direction between guide vanes
define(GVtc, 10)
// Number of cells in axial direction at guide vane inlet
define(GViac, 10)
// Number of cells in axial direction between guide vanes
define(GVbac, 10)
// Number of cells in axial direction at guide vane outlet
define(GVoac, 2)

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

// DRAFT TUBE REGION GEOMETRY AND MESH PROPERTIES
// "Draft tube" axial length (m)
define(DTal, 0.07)
// Number of sections per 360 degrees (integer!)
define(DTns, 5)
// Number of cells in radial direction in "draft tube"
define(DTrc, 10)
// Number of cells in tangential direction in "draft tube"
define(DTtc, 10)
// Number of cells in axial direction in "draft tube"
define(DTac, 7)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// TANGENTIAL PITCHES (RADIANS)
// Guide vane region
define(GVp, calc(2*pi/GVnb))
// Runner region
define(RUp, calc(2*pi/RUnb))
// Draft tube region
define(DTp, calc(2*pi/DTns))

// TANGENTIAL SHIFTS BETWEEN AXIAL LEVELS (BOTTOM-UP)
// Tangential shift from level DT0 to DT1
define(DTts01, calc(5*DTp))
// Runner region
// Tangential shift from level RU0 to RU1
define(RUts01, calc(-1/10*RUp))
// Tangential shift from level RU1 to RU2
define(RUts12, calc(-4/5*RUp))
// Tangential shift from level RU2 to RU3
define(RUts23, calc(-1/10*RUp))
// Guide vane region
// Tangential shift from level GV0 to GV1
define(GVts01, calc(1/10*GVp))
// Tangential shift from level GV1 to GV2
define(GVts12, calc(1/2*GVp))
// Tangential shift from level GV2 to GV3
define(GVts23, calc(0*GVp))

// AXIAL/TANGENTIAL BASE POINTS FOR EACH LEVEL (BOTTOM-UP):
// (CENTER OF RUNNER SET TO THETA=0, Z=0)
// Draft tube:
define(DTa0, calc(-RUoal-0.5*RUal-DTal)) //Center runner
define(DTt0, calc(-0.5*RUp-(0.5*RUts12)-(0*DTts01))) // Straight draft tube!
define(DTt1, calc(-0.5*RUp-(0.5*RUts12))) //Center runner
// Runner:
define(RUa0, calc(-RUoal-0.5*RUal)) //Center runner
define(RUt0, calc(-0.5*RUp-(0.5*RUts12))) //Center runner
define(RUt1, calc(RUt0+RUts01))
define(RUt2, calc(RUt1+RUts12))
define(RUt3, calc(RUt2+RUts23))
// Guide vane:
define(GVa0, calc(0.5*RUal+RUial)) //Center runner
define(GVt0, calc(-0.5*RUp-(0.5*RUts12)+RUts01+RUts12+RUts23)) //Center runner
define(GVt1, calc(GVt0+GVts01))
define(GVt2, calc(GVt1+GVts12))
define(GVt3, calc(GVt2+GVts23))

vertices //(radial [m], tangential [radians], axial [m])
(
//Guide vane hub:
    (hr GVt0 GVa0) vlabel(GV0lb)
    (hr calc(GVt0+GVp) GVa0) vlabel(GV0rb)
    (hr GVt1 calc(GVa0+GVoal)) vlabel(GV1lb)
    (hr calc(GVt1+GVp) calc(GVa0+GVoal)) vlabel(GV1rb)
    (hr GVt2 calc(GVa0+GVoal+GVbal)) vlabel(GV2lb)
    (hr calc(GVt2+GVp) calc(GVa0+GVoal+GVbal)) vlabel(GV2rb)
    (hr GVt3 calc(GVa0+GVoal+GVbal+GVial)) vlabel(GV3lb)
    (hr calc(GVt3+GVp) calc(GVa0+GVoal+GVbal+GVial)) vlabel(GV3rb)

//Guide vane shroud:
    (sr GVt0 GVa0) vlabel(GV0lt)
    (sr calc(GVt0+GVp) GVa0) vlabel(GV0rt)
    (sr GVt1 calc(GVa0+GVoal)) vlabel(GV1lt)
    (sr calc(GVt1+GVp) calc(GVa0+GVoal)) vlabel(GV1rt)
    (sr GVt2 calc(GVa0+GVoal+GVbal)) vlabel(GV2lt)
    (sr calc(GVt2+GVp) calc(GVa0+GVoal+GVbal)) vlabel(GV2rt)
    (sr GVt3 calc(GVa0+GVoal+GVbal+GVial)) vlabel(GV3lt)
    (sr calc(GVt3+GVp) calc(GVa0+GVoal+GVbal+GVial)) vlabel(GV3rt)

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

//Draft tube hub:
    (hr DTt0 DTa0) vlabel(DT0lb)
    (hr calc(DTt0+DTp) DTa0) vlabel(DT0rb)
    (hr DTt1 calc(DTa0+DTal)) vlabel(DT1lb)
    (hr calc(DTt1+DTp) calc(DTa0+DTal)) vlabel(DT1rb)

//Draft tube shroud:
    (sr DTt0 DTa0) vlabel(DT0lt)
    (sr calc(DTt0+DTp) DTa0) vlabel(DT0rt)
    (sr DTt1 calc(DTa0+DTal)) vlabel(DT1lt)
    (sr calc(DTt1+DTp) calc(DTa0+DTal)) vlabel(DT1rt)
);

blocks
(
//Guide vane:
    hex2D(GV0l, GV0r, GV1r, GV1l)
    (GVtc GVoac GVrc)
    simpleGrading (1 1 1)

    hex2D(GV1l, GV1r, GV2r, GV2l)
    (GVtc GVbac GVrc)
    simpleGrading (1 1 1)

    hex2D(GV2l, GV2r, GV3r, GV3l)
    (GVtc GViac GVrc)
    simpleGrading (1 1 1)

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

//Draft tube:
    hex2D(DT0l, DT0r, DT1r, DT1l)
    (DTtc DTac DTrc)
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
//Guide vane
    spline GV1lt GV2lt
    (
        (sr calc(GVt1+0.75*(GVt2-(GVt1))) calc(GVa0+GVoal+0.5*GVbal))
    )
    spline GV1lb GV2lb
    (
        (hr calc(GVt1+0.75*(GVt2-(GVt1))) calc(GVa0+GVoal+0.5*GVbal))
    )
    spline GV1rt GV2rt
    (
        (sr calc(GVt1+GVp+0.65*(GVt2-(GVt1))) calc(GVa0+GVoal+0.5*GVbal))
    )
    spline GV1rb GV2rb
    (
        (hr calc(GVt1+GVp+0.65*(GVt2-(GVt1))) calc(GVa0+GVoal+0.5*GVbal))
    )
);

boundary
(
    GVINLET
    {
        type            patch;
        faces
        (
            quad2D(GV3r, GV3l)
        );
    }

    GVOUTLET
    {
        type            overlapGgi;
        shadowPatch     RUINLET;
        zone            GVOUTLETZone;
        rotationAxis    ( 0 0 1 );
        nCopies         5;
        faces
        (
            quad2D(GV0l, GV0r)
        );
    }

    GVCYCLIC1
    {
        type             cyclicGgi;
        shadowPatch      GVCYCLIC2;
        zone             GVCYCLIC1Zone;
        bridgeOverlap    false;
        rotationAxis     (0 0 1);
        rotationAngle    72;
        separationOffset (0 0 0);
        faces
        (
            quad2D(GV1l, GV0l)
            quad2D(GV3l, GV2l)
        );
    }

    GVCYCLIC2
    {
        type             cyclicGgi;
        shadowPatch      GVCYCLIC1;
        zone             GVCYCLIC2Zone;
        bridgeOverlap    false;
        rotationAxis     (0 0 1);
        rotationAngle    -72;
        separationOffset (0 0 0);
        faces
        (
            quad2D(GV0r, GV1r)
            quad2D(GV2r, GV3r)
        );
    }

    //GVCYCLIC
    //{
    //    type cyclic;
    //    faces
    //    (
    //        quad2D(GV1l, GV0l)
    //        quad2D(GV3l, GV2l)
    //        quad2D(GV0r, GV1r)
    //        quad2D(GV2r, GV3r)
    //    );
    //}

    GVBLADE
    {
        type wall;
        faces
        (
            quad2D(GV2l, GV1l)
            quad2D(GV1r, GV2r)
        );
    }

    GVHUB
    {
        type wall;
        faces
        (
            backQuad(GV0l, GV0r, GV1r, GV1l)
            backQuad(GV1l, GV1r, GV2r, GV2l)
            backQuad(GV2l, GV2r, GV3r, GV3l)
        );
    }

    GVSHROUD
    {
        type wall;
        faces
        (
            frontQuad(GV0l, GV0r, GV1r, GV1l)
            frontQuad(GV1l, GV1r, GV2r, GV2l)
            frontQuad(GV2l, GV2r, GV3r, GV3l)
        );
    }

    RUINLET
    {
        type            overlapGgi;
        shadowPatch     GVOUTLET;
        zone            RUINLETZone;
        rotationAxis    ( 0 0 1 );
        nCopies         5;
        faces
        (
            quad2D(RU3r, RU3l)
        );
    }

    RUOUTLET
    {
        type            overlapGgi;
        shadowPatch     DTINLET;
        zone            RUOUTLETZone;
        rotationAxis    ( 0 0 1 );
        nCopies         5;
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

    //RUCYCLIC
    //{
    //    type cyclic;
    //    faces
    //    (
    //        quad2D(RU1l, RU0l)
    //        quad2D(RU3l, RU2l)
    //        quad2D(RU0r, RU1r)
    //        quad2D(RU2r, RU3r)
    //    );
    //}

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

    DTINLET
    {
        type            overlapGgi;
        shadowPatch     RUOUTLET;
        zone            DTINLETZone;
        rotationAxis    ( 0 0 1 );
        nCopies         5;
        faces
        (
            quad2D(DT1r, DT1l)
        );
    }

    DTOUTLET
    {
        type patch;
        faces
        (
            quad2D(DT0l, DT0r)
        );
    }

    DTCYCLIC1
    {
        type             cyclicGgi;
        shadowPatch      DTCYCLIC2;
        zone             DTCYCLIC1Zone;
        bridgeOverlap    false;
        rotationAxis     (0 0 1);
        rotationAngle    72;
        separationOffset (0 0 0);
        faces
        (
            quad2D(DT1l, DT0l)
        );
    }

    DTCYCLIC2
    {
        type             cyclicGgi;
        shadowPatch      DTCYCLIC1;
        zone             DTCYCLIC2Zone;
        bridgeOverlap    false;
        rotationAxis     (0 0 1);
        rotationAngle    -72;
        separationOffset (0 0 0);
        faces
        (
            quad2D(DT0r, DT1r)
        );
    }

    //DTCYCLIC
    //{
    //    type cyclic;
    //    faces
    //    (
    //        quad2D(DT1l, DT0l)
    //        quad2D(DT0r, DT1r)
    //    );
    //}

    DTHUB
    {
        type wall;
        faces
        (
            backQuad(DT0l, DT0r, DT1r, DT1l)
        );
    }

    DTSHROUD
    {
        type wall;
        faces
        (
            frontQuad(DT0l, DT0r, DT1r, DT1l)
        );
    }
);


// ************************************************************************* //
