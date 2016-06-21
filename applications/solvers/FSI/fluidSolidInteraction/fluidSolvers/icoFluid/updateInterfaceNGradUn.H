{
    Info << "Update normal derivative of normal velocity" << endl;

    volVectorField phiU = fvc::reconstruct(phi_);

    vectorField n = mesh.boundary()[fsi.fluidPatchIndex()].nf();

    scalarField UnP =
        (n&phiU.boundaryField()[fsi.fluidPatchIndex()].patchInternalField());

    scalarField UnInt =
        phi_.boundaryField()[fsi.fluidPatchIndex()]
       /mesh.magSf().boundaryField()[fsi.fluidPatchIndex()];

    nGradUn_ =
        (UnInt - UnP)
       *mesh.deltaCoeffs().boundaryField()[fsi.fluidPatchIndex()];

    Info << "nGradUn, max: " << gMax(nGradUn_) << ", min: "
        << gMin(nGradUn_) << ", avg: " << gAverage(nGradUn_) << endl;

//     bool secondOrderCorrection = true;

//     if (secondOrderCorrection)
//     {
//         // Correct normal component of phiU
//         // befor gradient calculation
//         forAll(phiU.boundaryField(), patchI)
//         {
//                 vectorField n =
//               mesh().Sf().boundaryField()[patchI]
//               /mesh().magSf().boundaryField()[patchI];

//             phiU.boundaryField()[patchI] +=
//                n
//               *(
//                    (
//                        phi_.boundaryField()[patchI]
//                       /mesh().magSf().boundaryField()[patchI]
//                    )
//                  - (n&phiU.boundaryField()[patchI])
//                );
//         }

//         // Calc gradient
//         tensorField gradPhiUp =
//             fvc::grad(phiU)().boundaryField()[aPatchID()].patchInternalField();
//             nGradUn() = 2*nGradUn() - (nA&(gradPhiUp&nA));
//     }
}
