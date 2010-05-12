scalar velCmptMag = sqrt(moleculeCloud::kb*temperature/mass);

for (molN = totalMols; molN < totalMols + totalZoneMols; molN++)
{
    // Assign velocity: random direction, magnitude determined by desired
    // maxwellian distribution at temperature

    // Temperature gradients could be created by specifying a gradient in the
    // zone subDict, or by reading a field from a mesh.

    // The velocities are treated on a zone-by-zone basis for the purposes of
    // removal of bulk momentum - hence nMols becomes totalZoneMols

    velocity = vector
    (
        velCmptMag*rand.GaussNormal(),
        velCmptMag*rand.GaussNormal(),
        velCmptMag*rand.GaussNormal()
    );

    momentumSum += mass*velocity;

    initialVelocities.append(velocity);
}


