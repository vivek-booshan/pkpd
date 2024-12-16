function dydt = ihatemylifeODE(t, y, p)
    
    GI_chol = y(1);
    peripheral_chol = y(2);
    liver_chol = y(3);
    clearance_chol = y(4);

    GI_TG = y(5);
    peripheral_TG = y(6);
    liver_TG = y(7);
    clearance_TG = y(8);

    ROS = y(9);
    oxLDL = y(10);

    oral_statin = y(11);
    statin = y(12);
    liver_statin = y(13);
    clearance_statin = y(14);

    dGI_chol = -p.V.GI * p.GI.chylomicron.chol * GI_chol;
    dGI_TG = -p.V.GI * p.GI.chylomicron.TG * GI_TG;

    dperipheral_chol = (...
        + p.V.GI * p.GI.chylomicron.chol * GI_chol ...
        + p.V.liver * p.liver.LDL.chol * liver_chol * log(ROS) ...
        + p.V.liver * p.liver.HDL.chol * liver_chol ... 
        - p.V.peripheral * p.peripheral.LDL.chol * peripheral_chol * log(liver_statin) ...
        - p.V.peripheral * p.peripheral.chylomicron.chol * peripheral_chol ...
        - p.V.peripheral * p.peripheral.HDL.chol * peripheral_chol * (1/ROS) ...
        - p.V.peripheral * p.peripheral.clearance.chol * peripheral_chol ...
        ) / p.V.peripheral;
    dperipheral_TG = (...
        + p.V.GI * p.GI.chylomicron.TG * GI_TG + ...
        + p.V.liver * p.liver.LDL.chol * liver_chol ...
        + p.V.liver * p.liver.HDL.chol * liver_chol ...
        - p.V.peripheral * p.peripheral.chylomicron.TG * peripheral_TG ...
        - p.V.peripheral * p.peripheral.HDL.TG * peripheral_TG ...
        - p.V.peripheral * p.peripheral.clearance.TG * peripheral_TG...
        ) / p.V.peripheral;
    
    dliver_chol = (...
        + p.V.peripheral * p.peripheral.chylomicron.chol * peripheral_chol ...
        + p.V.peripheral * p.peripheral.HDL.chol * peripheral_chol * (1/ROS) ...
        + p.V.peripheral * p.peripheral.clearance.chol * peripheral_chol ...
        + p.V.peripheral * p.peripheral.LDL.chol * peripheral_chol * log(liver_statin) ...
        - p.V.liver * p.liver.LDL.chol * liver_chol * log(ROS) ...
        - p.V.liver * p.liver.HDL.chol * liver_chol ...
        - p.V.liver * p.liver.clearance.chol * liver_chol ...
        - p.V.liver * p.liver.oxLDL.chol * liver_chol * log(ROS) ...
        ) / p.V.liver;

    dliver_TG = (...
        + p.V.peripheral * p.peripheral.chylomicron.TG * peripheral_TG ...
        + p.V.peripheral * p.peripheral.HDL.TG * peripheral_TG ...
        + p.V.peripheral * p.peripheral.clearance.TG * peripheral_TG ...
        - p.V.liver * p.liver.LDL.TG * liver_TG ...
        - p.V.liver * p.liver.HDL.TG * liver_TG ...
        - p.V.liver * p.liver.clearance.TG * liver_TG ...
        ) / p.V.liver;

    dclearance_chol = (...
        + p.V.peripheral * p.peripheral.clearance.chol * peripheral_chol ...
        + p.V.liver * p.liver.clearance.chol * liver_chol ...
        );
    dclearance_TG = (...
        + p.V.peripheral * p.peripheral.clearance.TG * peripheral_TG ...
        + p.V.liver * p.liver.clearance.TG * liver_TG ...
        );
    dROS = (...
        + p.basalROS * ROS * p.foodProductionMultiplier ...
        - p.antioxidant * ROS * p.foodAntiOxMultiplier ...
        );

    doxLDL = p.V.liver * p.liver.oxLDL.chol * liver_chol * ROS; 
    
    doral_statin = - p.F * p.ka * oral_statin;
    dstatin = + p.F * p.ka * oral_statin - p.statin_liver * statin;
    dliver_statin = + p.liver_statin * statin - p.liver.clearance.statin * liver_statin;
    dclearance_statin = + p.liver.clearance.statin * liver_statin;

    dydt = [...
        dGI_chol; dperipheral_chol; dliver_chol; dclearance_chol; ...
        dGI_TG; dperipheral_TG; dliver_TG; dclearance_TG; ...
        dROS; doxLDL ...
        ];
end