    forAll (cr, crI)
    {
        UEqn.setComponentReference
        (
            cr[crI].patchIndex(),
            cr[crI].faceIndex(),
            cr[crI].dir(),
            cr[crI].value()
        );
    }
