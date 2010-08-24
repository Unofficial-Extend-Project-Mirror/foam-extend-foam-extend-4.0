{
    Info << "moving piston points" << endl;

    boolList scaleDisp(nPoints(), true);
    label nScaled = nPoints();

    List<bool> pistonPoint(newPoints.size(), false);

    label pistonPtsIndex = pointZones().findZoneID("pistonPoints");
    const labelList& pistonPoints = pointZones()[pistonPtsIndex];

    const scalarField& movingPointsMPiston = movingPointsMaskPiston();

    forAll(pistonPoints, i)
    {
        label pointI = pistonPoints[i];
        pistonPoint[pointI] = true;
        point& p = newPoints[pointI];
    }

    // Always move piston
    scalar pistonTopZ = -GREAT;

    forAll(pistonPoints, i)
    {
        point& p = newPoints[pistonPoints[i]];
        p.z() = p.z() + deltaZ;
    }
}
