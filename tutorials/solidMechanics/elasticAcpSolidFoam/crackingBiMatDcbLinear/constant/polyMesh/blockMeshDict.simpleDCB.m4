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
    location    "constant/polyMesh";
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Model Description
// Plane strain 2-D model of simple DCB
// two beams adhered together
// thickness of each beam can be set
// or the material of each beam could be set using setFields

// Setup m4 stuff
changecom(//)changequote([,])
define(calc, [esyscmd(perl -e 'printf ($1)')])
//define(VCOUNT, 0)
//define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])

// define geometry in mm
define(bl, 60) // beams length
define(bh, 10) // individual beam height
define(pn, 10) // pre-notch length
define(z, 1) // out of plane thickness

// define mesh density
define(meshPn, 10) // number of cells in beam length direction before notch
define(meshBl, 60) // number of cells in beam length direction after notch
define(meshBh, 10)  // number of cells in beam height direction

// start of blockMeshDict

convertToMeters 0.001;

vertices
(
 //- dimension in mm
 (0 -bh 0)
 (pn -bh 0)
 (bl -bh 0)
 (bl 0 0)
 (bl bh 0)
 (pn bh 0)
 (0 bh 0)
 (0 0 0)
 (pn 0 0)
 (0 0 0) // 9

 (0 -bh z)
 (pn -bh z)
 (bl -bh z)
 (bl 0 z)
 (bl bh z)
 (pn bh z)
 (0 bh z)
 (0 0 z)
 (pn 0 z)
 (0 0 z) // 19

);

blocks
(
//  hex (0 1 8 9 10 11 18 19) sheet (meshPn meshBh 1) simpleGrading (1 1 1)
//  hex (1 2 3 8 11 12 13 18) sheet (meshBl meshBh 1) simpleGrading (1 1 1)
//  hex (7 8 5 6 17 18 15 16) sheet (meshPn meshBh 1) simpleGrading (1 1 1)
//  hex (8 3 4 5 18 13 14 15) sheet (meshBl meshBh 1) simpleGrading (1 1 1)
 hex (0 1 8 9 10 11 18 19) (meshPn meshBh 1) simpleGrading (1 1 1)
 hex (1 2 3 8 11 12 13 18) (meshBl meshBh 1) simpleGrading (1 1 1)
 hex (7 8 5 6 17 18 15 16) (meshPn meshBh 1) simpleGrading (1 1 1)
 hex (8 3 4 5 18 13 14 15) (meshBl meshBh 1) simpleGrading (1 1 1)
 );

edges
(
);

patches
(
 patch tractionFree
 (
  (15 16 6 5)
  (14 15 5 4)

  (0 1 11 10)
  (1 2 12 11)

  (2 3 13 12)
  (3 4 14 13)

  (7 8 18 17)
  (19 18 8 9)
  )

 patch topLoading
 (
  (17 16 6 7)
  )

 patch bottomLoading
 (
  (10 19 9 0)
  )

 empty back
 (
  (7 6 5 8)
  (8 5 4 3)
  (0 9 8 1)
  (1 8 3 2)
  )

 empty front
 (
  (10 11 18 19)
  (11 12 13 18)
  (17 18 15 16)
  (18 13 14 15)
  )

 cohesive crack ()
 );

mergePatchPairs
(
);

// ************************************************************************* //
