scalar initVelMag =
    sqrt
    (
        3.0*(1.0 - 1.0 / totalZoneMols)
       *moleculeCloud::kb*temperature
       /mass
    );

for (molN = totalMols; molN < totalMols + totalZoneMols; molN++)
{
    // Assign velocity: random direction, magnitude determined by desired
    // temperature

    // Temperature gradients could be created by specifying a gradient in the
    // zone subDict, or by reading a field from a mesh.

    // The velocities are treated on a zone-by-zone basis for the purposes of
    // removal of bulk momentum - hence nMols becomes totalZoneMols

    velocity = (2.0*rand.vector01() - vector::one);

    velocity *= initVelMag/mag(velocity);

    momentumSum += mass*velocity;

    initialVelocities.append(velocity);
}
