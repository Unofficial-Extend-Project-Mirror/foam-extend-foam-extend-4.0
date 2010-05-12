// use tetrahedral decomposition of the engine mesh

{
    tetPointVectorField& motionU = mSolver.motionU();

    vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

    Info << " Set valve velocity " << endl;
    
    forAll (valves_, valveI)
    {
        // If valve is present in geometry, set the motion
       
       Info << " Set valve Bottom velocity " << endl;
       
       if (valves_[valveI].bottomPatchID().active())
        {
            vector valveVel =
                -valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

            // Bottom of the valve moves with given velocity
            motionU.boundaryField()[valves_[valveI].bottomPatchID().index()] ==
                valveVel;

            if (debug)
            {
                Info<< "Valve " << valveI << " lift: "
                    << valves_[valveI].curLift()
                    << " velocity: " << valves_[valveI].curVelocity()
                    << endl;
            }
        }

       
       Info << " Set valve poppet velocity " << endl;
       
        if (valves_[valveI].poppetPatchID().active())
        {

            // Top of the valve does not move
            motionU.boundaryField()[valves_[valveI].poppetPatchID().index()] ==
                valveVel;
        }

       
       Info << " Set valve curInCyl velocity " << endl;
       
        if (valves_[valveI].curtainInCylinderPatchID().active())
        {

            label cicIndex = valves_[valveI].curtainInCylinderPatchID().index();
            
            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[cicIndex]
                );
            
            scalarField curtainInCylinderPoints =
                            motionU.boundaryField()[cicIndex].patch()
                            .localPoints().component(vector::Z);

            pp.refValue() = pistonVel;
            pp.valueFraction() = 
            (1.0-pos(curtainInCylinderPoints - pistonPosition()))*vector(0,0,1);
                
            
                            
        }
       
       Info << " Set valve curInport velocity " << endl;
       

        if (valves_[valveI].curtainInPortPatchID().active())
        {
            label cipIndex = valves_[valveI].curtainInPortPatchID().index();
            
            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[cipIndex]
                );
            
            scalarField curtainInPortPoints =
                            motionU.boundaryField()[cipIndex].patch()
                            .localPoints().component(vector::Z);

            pp.refValue() = pistonVel;
            pp.valueFraction() = 
            (1.0-pos(curtainInPortPoints - pistonPosition()))*vector(0,0,1);;

        }
        
       
       Info << " Set valve stem velocity " << endl;
       
        if (valves_[valveI].stemPatchID().active())
        {
            label stemIndex = valves_[valveI].stemPatchID().index();
            
            componentMixedTetPolyPatchVectorField& pp =
                refCast<componentMixedTetPolyPatchVectorField>
                (
                    motionU.boundaryField()[stemIndex]
                );
            
            pp.refValue() = vector::zero;
            pp.valueFraction() = vector(0,0,1);
        }
    }    
/*        
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

            pp.refValue() = pistonVel;
            pp.valueFraction() = 
            (1.0-pos(linerPoints - pistonPosition()))*vector(0,0,1);;

        }
        
    

    // Set piston velocity
    if (piston().patchID().active())
    {

        if (debug)
        {
            Info<< "Piston velocity: " << pistonVel;
        }

        componentMixedTetPolyPatchVectorField& pp =
            refCast<componentMixedTetPolyPatchVectorField>
            (
                motionU.boundaryField()[piston().patchID().index()]
            );

        pp.refValue() = pistonVel;
        
    }
*/    
    
    
    mSolver.solve();
	
    Info << "msolved" << endl;
    
//  Setting the boundary position

    {
        label cylinderHeadIndex = boundaryMesh().findPatchID("cylinderHead");
        
        // Top of the valve does not move
        motionU.boundaryField()[cylinderHeadIndex] ==
            vector::zero;
    }
        
    
}
