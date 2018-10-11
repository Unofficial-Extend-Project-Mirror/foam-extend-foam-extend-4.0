{
    // Make the fluxes absolute (using the ddt(U) scheme)
    fvc::makeAbsolute(phi, U);

    // Motion inside PIMPLE loop: simple IB update
    bool meshChanged = mesh.update();

#   include "updateIbPatchFields.H"
#   include "updateIbMasks.H"
#   include "volContinuity.H"

    if (runTime.outputTime())
    {
        volScalarField divMeshPhi
        (
            "divMeshPhi",
            -fvc::div(mesh.phi())
        );

        divMeshPhi.internalField() +=
            (1.0 - mesh.V0()/mesh.V())/runTime.deltaT().value();

        divMeshPhi = mag(divMeshPhi);
        divMeshPhi.write();

        U.write();
        mesh.phi().write();
    }

    // Fluxes will be corrected to absolute velocity
    // HJ, 6/Feb/2009
#   include "correctPhi.H"

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

    if (checkMeshCourantNo)
    {
#       include "immersedBoundaryCourantNo.H"
    }

    // Correct velocity boundary conditions after mesh motion, with relative
    // fluxes
    U.correctBoundaryConditions();
}
