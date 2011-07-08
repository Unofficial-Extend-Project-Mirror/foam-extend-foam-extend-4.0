
    vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

{
    tetPointVectorField& motionU = mSolver.motionU();

    Info << "setting the piston velocity" << endl;

    if (piston().patchID().active())
    {

        if (debug)
        {
            Info << "Piston velocity: " << pistonVel;
        }

        motionU.boundaryField()[piston().patchID().index()] == pistonVel;

    }

}
