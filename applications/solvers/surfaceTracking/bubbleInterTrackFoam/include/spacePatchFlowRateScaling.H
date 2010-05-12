if(spacePatchID != -1)
{
    fvsPatchScalarField& spacePhi = phi.boundaryField()[spacePatchID];

    scalar outflowScaling =
    (
      - sum(neg(spacePhi)*spacePhi) + SMALL
    )/
    (
        sum(pos(spacePhi)*spacePhi) + SMALL
    );

    spacePhi += pos(spacePhi)*spacePhi*(outflowScaling - 1.0);

    U.boundaryField()[spacePatchID] +=
        pos(spacePhi)*U.boundaryField()[spacePatchID]*(outflowScaling - 1.0);
}
