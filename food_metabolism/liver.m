% TODO : minor optimizations involving passing by reference to decrease
% overhead from dydt initialization and hopefully a static inline to reduce
% function overhead (need to figure out how to compile this a c code)


% TODO : insulin needs to be self regulating
% TODO : glucagon needs to be self regulating
% TODO : somatostatin needs to be self regulating
% TODO : redefine glucose to modulate insulin, glucagon, and somatostatin secretion
% TODO : need to add clearances
% 
function dydt = liver(t, y, p)
    dydt = zeros(4, 1);
    dinsulin = insulin(t, y, p);
    dglucagon = glucagon(t, y, p);
    dsomatostatin = somatostatin(t, y, p);
    dglucose = glucose(t, y, p);

    dydt = dinsulin + dglucagon + dsomatostatin + dglucose;
end

function dydt = insulin(t, y, p)
    dydt = zeros(4, 1);
    dydt(1) = (...
        + p.insulin_pancreas_portal * y(1) * p.V_pancreas ...
        
        + p.insulin_blood_liver * y(1) * p.V_blood ... 
        - p.insulin_liver_blood * y(1) * p.V_liver ...
    ) / p.V_liver;
end

function dydt = glucagon(t, y, p)
    dydt = zeros(4, 1);
    dydt(2) = (...
        + p.glucagon_pancreas_portal * y(2) * p.V_pancreas ...
        
        + p.glucagon_blood_liver * y(2) * p.V_blood ... 
        - p.glucagon_liver_blood * y(2) * p.V_liver ...
    ) / p.V_liver;
end

function dydt = somatostatin(t, y, p)
    dydt = zeros(4, 1); 
    dydt(3) = (...
        + p.somatostatin_pancreas_portal * y(3) * p.V_pancreas ...
        
        + p.somatostatin_blood_liver * y(3) * p.V_blood ... 
        - p.somatostatin_liver_blood * y(3) * p.V_liver ...
    ) / p.V_liver;
end

function dydt = glucose(t, y, p)
    dydt = zeros(4, 1);
    dydt(4) = (...
        + p.glucose_pancreas_portal * y(4) * p.V_pancreas ...
        
        + p.glucose_blood_liver * y(4) * p.V_blood ... 
        - p.glucose_liver_blood * y(4) * p.V_liver ...
    ) / p.V_liver;
end

