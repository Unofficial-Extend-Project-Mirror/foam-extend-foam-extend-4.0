if(runTime.timeIndex() > 1)
{
    fvsPatchScalarField& spacePhi = phi.boundaryField()[spacePatchID];

    scalar inletFlux = gSum(neg(spacePhi)*spacePhi);

    scalar outletFlux = gSum(pos(spacePhi)*spacePhi);

    if(outletFlux < VSMALL)
    {
        outletFlux = VSMALL;
    }

    scalar outflowScaling = -inletFlux/outletFlux;

    spacePhi += pos(spacePhi)*spacePhi*(outflowScaling - 1.0);

//     Pout << "Space patch continuity: " 
//         << gSum(phi.boundaryField()[spacePatchID]) << endl;

//     U.boundaryField()[spacePatchID] ==
//         U.boundaryField()[spacePatchID]
//       + pos(spacePhi)*U.boundaryField()[spacePatchID]*(outflowScaling - 1.0);
}
