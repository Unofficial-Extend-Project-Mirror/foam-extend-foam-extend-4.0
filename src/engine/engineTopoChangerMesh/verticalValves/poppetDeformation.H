if(poppetDeformation)
{

    Info << "Poppet valve deformation" << endl;

    tetPointVectorField& motionU = mSolver.motionU();
    // Blocking vertical motion

    motionU.internalField().replace(vector::Z, 0);
    
    // Set valve velocity
    forAll (valves_, valveI)
    {

        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();
        
        scalar zPoppetMin = -GREAT;
        scalar zPoppetMax = GREAT;
        scalar zCurtainMax = GREAT;
                
        // If valve is present in geometry, set the motion
        if (valves_[valveI].bottomPatchID().active())
        {
            // Bottom of the valve moves with given velocity
            motionU.boundaryField()[valves_[valveI].bottomPatchID().index()] ==
                vector::zero;
        }
        
        if (valves_[valveI].poppetPatchID().active()) 
        {
            
            zPoppetMin = min
            (
                motionU.boundaryField()
                [valves_[valveI].poppetPatchID().index()].patch()
                .localPoints().component(vector::Z)
            );
            zPoppetMax = max
            (
                motionU.boundaryField()
                [valves_[valveI].poppetPatchID().index()].patch()
                .localPoints().component(vector::Z)
            );
            
            if (valves_[valveI].curLift() > valves_[valveI].deformationLift())
            {
                // Top of the valve does not move
                motionU.boundaryField()
                    [valves_[valveI].poppetPatchID().index()] ==
                    vector::zero;
            }
            else
            {
                motionU.boundaryField()
                    [valves_[valveI].poppetPatchID().index()] ==
                    valveVel;
            }
        
        }

        if (valves_[valveI].curtainInCylinderPatchID().active())
        {

            label cicIndex = valves_[valveI].curtainInCylinderPatchID().index();
            
            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[cicIndex]
                );
            
            pp.refValue() = vector::zero;
            pp.valueFraction() = vector::zero;
                            
        }
        
        scalar zDetMin = -GREAT;        
        scalar zDetMax = GREAT;        

        if (valves_[valveI].detachInPortPatchID().active())
        {            
            zDetMin = min
            (
                motionU.boundaryField()
                [valves_[valveI].detachInPortPatchID().index()].patch()
                .localPoints().component(vector::Z)
            );
            
            zDetMax = max
            (
                motionU.boundaryField()
                [valves_[valveI].detachInPortPatchID().index()].patch()
                .localPoints().component(vector::Z)
            );
        
            motionU.boundaryField()
                [valves_[valveI].detachInPortPatchID().index()] ==
                vector::zero;            
        }

        if (valves_[valveI].detachInCylinderPatchID().active())
        {
            motionU.boundaryField()
                [valves_[valveI].detachInCylinderPatchID().index()] ==
                vector::zero;
        }
        

        if (valves_[valveI].curtainInPortPatchID().active())
        {
            label cipIndex = valves_[valveI].curtainInPortPatchID().index();

            zCurtainMax = max
            (
                motionU.boundaryField()[cipIndex].patch()
                .localPoints().component(vector::Z)
            );
            
            if (valves_[valveI].curLift() > valves_[valveI].deformationLift())
            {
                componentMixedTetPolyPatchVectorField& pp =
                    refCast<componentMixedTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[cipIndex]
                    );
            
                pp.refValue() = vector::zero;
                pp.valueFraction() = vector::zero;
            }
            else
            {
            
//              poppet deformation
            
                componentMixedTetPolyPatchVectorField& pp =
                    refCast<componentMixedTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[cipIndex]
                    );
                
                scalarField cicPoints = 
                motionU.boundaryField()[cipIndex].patch()
                .localPoints().component(vector::Z);
                
                pp.valueFraction() = pos(cicPoints - zPoppetMin)*vector(0,0,1);

                pp.refValue() =
                    (cicPoints - zCurtainMax)/
                    (zPoppetMin - zCurtainMax)*valveVel;            
            }

        }

        if (valves_[valveI].stemPatchID().active())
        {
            label stemIndex = valves_[valveI].stemPatchID().index();
            
            if (valves_[valveI].curLift() > valves_[valveI].deformationLift())
            {
                componentMixedTetPolyPatchVectorField& pp =
                    refCast<componentMixedTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[stemIndex]
                    );
            
                pp.refValue() = vector::zero;
                pp.valueFraction() = vector::zero;
            }
            else
            {
            
//              poppet deformation
            
                componentMixedTetPolyPatchVectorField& pp =
                    refCast<componentMixedTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[stemIndex]
                    );
                
                scalarField stemPoints = 
                motionU.boundaryField()[stemIndex].patch()
                .localPoints().component(vector::Z);

                pp.valueFraction() =
                    (1.0 - pos(stemPoints - zPoppetMax))*vector(0,0,1);

                pp.refValue() =
                    (stemPoints - zDetMax)/(zPoppetMax - zDetMax)*valveVel;
            }
        }


        {
            label linerIndex = boundaryMesh().findPatchID("liner");
            
            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[linerIndex]
                );
            
            scalarField linerPoints =
                            motionU.boundaryField()[linerIndex].patch()
                            .localPoints().component(vector::Z);

            pp.refValue() = vector::zero;
            pp.valueFraction() = vector::zero;

        }
        
    }

    // Set piston velocity
    if (piston().patchID().active())
    {

        componentMixedTetPolyPatchVectorField& pp =
            refCast<componentMixedTetPolyPatchVectorField>
            (
                motionU.boundaryField()[piston().patchID().index()]
            );

        pp.refValue() = vector::zero;
        
    }

    {
        label cylinderHeadIndex = boundaryMesh().findPatchID("cylinderHead");
        
        // Top of the valve does not move
        motionU.boundaryField()[cylinderHeadIndex] ==
            vector::zero;
    }
    
    mSolver.solve();
    
//  Setting the boundary position
        
    movePoints(mSolver.curPoints());

    // Blocking vertical motion
    motionU.internalField().replace(vector::Z, 0);
}
