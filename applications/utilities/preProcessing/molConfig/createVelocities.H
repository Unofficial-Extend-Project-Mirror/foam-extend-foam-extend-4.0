vector velocity;

vector momentumSum = vector::zero;

if (velocityDistribution == "uniform")
{
#    include "uniform.H"
}

if (velocityDistribution == "maxwellian")
{
#    include "maxwellian.H"
}
