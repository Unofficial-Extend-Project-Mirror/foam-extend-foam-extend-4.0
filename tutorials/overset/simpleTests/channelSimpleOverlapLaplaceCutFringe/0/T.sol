/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     3.1                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  | For copyright notice see file Copyright         |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "1";
    object      T;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 1 0 0];

internalField   nonuniform List<scalar> 
60
(
1.05
1.15
1.25
1.35
1.45
1.55
1.05
1.15
1.25
1.35
1.45
1.55
1.05
1.15
1.25
1.35
1.45
1.55
1.05
1.15
1.25
1.35
1.45
1.55
1.05
1.15
1.25
1.35
1.45
1.55
1.45
1.55
1.65
1.75
1.85
1.95
1.45
1.55
1.65
1.75
1.85
1.95
1.45
1.55
1.65
1.75
1.85
1.95
1.45
1.55
1.65
1.75
1.85
1.95
1.45
1.55
1.65
1.75
1.85
1.95
)
;

boundaryField
{
    oversetFaces
    {
        type            overset;
        coupledFringe   yes;
        setHoleCellValue yes;
        holeCellValue   0;
        value           nonuniform 0();
    }
    left
    {
        type            fixedValue;
        value           uniform 1;
    }
    leftEnd
    {
        type            zeroGradient;
    }
    rightStart
    {
        type            zeroGradient;
    }
    right
    {
        type            fixedValue;
        value           uniform 2;
    }
    top
    {
        type            zeroGradient;
    }
    bottom
    {
        type            zeroGradient;
    }
    frontAndBack
    {
        type            empty;
    }
}


// ************************************************************************* //
