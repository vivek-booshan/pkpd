function dydt = pcsk9ODE(t, y, p)
    % State variables
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

    % PCSK9 inhibitor effect: Increased LDL clearance rate
    pcsk9_LDL_clearance = 1.5; % Increases LDL clearance by 50%

    % Differential equations
    dGI_chol = -p.V.GI * p.GI.chylomicron.chol * GI_chol;
    dGI_TG = -p.V.GI * p.GI.chylomicron.TG * GI_TG;

    dperipheral_chol = (p.V.GI * p.GI.chylomicron.chol * GI_chol ...
        - p.V.peripheral * p.peripheral.LDL.chol * peripheral_chol * pcsk9_LDL_clearance ...
        - p.V.peripheral * p.peripheral.clearance.chol * peripheral_chol) / p.V.peripheral;

    dperipheral_TG = (p.V.GI * p.GI.chylomicron.TG * GI_TG ...
        - p.V.peripheral * p.peripheral.clearance.TG * peripheral_TG) / p.V.peripheral;

    dliver_chol = (p.V.peripheral * p.peripheral.LDL.chol * peripheral_chol * pcsk9_LDL_clearance ...
        - p.V.liver * p.liver.clearance.chol * liver_chol) / p.V.liver;

    dliver_TG = (p.V.peripheral * p.peripheral.clearance.TG * peripheral_TG ...
        - p.V.liver * p.liver.clearance.TG * liver_TG) / p.V.liver;

    dclearance_chol = (p.V.peripheral * p.peripheral.clearance.chol * peripheral_chol ...
        + p.V.liver * p.liver.clearance.chol * liver_chol);

    dclearance_TG = (p.V.peripheral * p.peripheral.clearance.TG * peripheral_TG ...
        + p.V.liver * p.liver.clearance.TG * liver_TG);

    dROS = p.basalROS * ROS * p.foodProductionMultiplier ...
         - p.antioxidant * ROS * p.foodAntiOxMultiplier;

    doxLDL = p.V.liver * p.liver.oxLDL.chol * liver_chol * ROS;

    % Combine outputs
    dydt = [dGI_chol; dperipheral_chol; dliver_chol; dclearance_chol; ...
            dGI_TG; dperipheral_TG; dliver_TG; dclearance_TG; ...
            dROS; doxLDL];
end
