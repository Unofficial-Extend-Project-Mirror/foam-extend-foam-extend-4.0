{
  scalar maxDU = gMax(mag(DU.internalField()));

    relativeResidual =
        gMax
      (
       mag(DU.internalField() - DU.prevIter().internalField())
       /(maxDU + SMALL)
       );

    resid = (DU - DU.prevIter())/(mag(DU)+dimensionedScalar("small",dimLength,SMALL));
}

// {
//     scalarField magDU = mag(DU.internalField());

//     forAll(magDU, cellI)
//     {
//         if (magDU[cellI] < SMALL)
//         {
//             magDU[cellI] = SMALL;
//         }
//     }

//     relativeResidual =
//         gMax
//         (
//             mag
//             (
//                 DU.internalField()
//               - DU.prevIter().internalField()
//             )
//            /magDU
//         );
// }
