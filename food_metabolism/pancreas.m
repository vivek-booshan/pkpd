% TODO : minor optimizations involving passing by reference to decrease
% overhead from dydt initialization and hopefully a static inline to reduce
% function overhead (need to figure out how to compile this a c code)

% TODO : insulin needs to be self regulating
% TODO : glucagon needs to be self regulating
% TODO : somatostatin needs to be self regulating
% TODO : redefine glucose to modulate insulin, glucagon, and somatostatin secretion

function dydt = pancreas(t, y, p)
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
        + p.insulin_blood_pancreas * y(1) * p.V_blood ...
        - p.insulin_pancreas_portal * y(1) * p.V_pancreas ...
    ) / p.V_pancreas;
end

function dydt = glucagon(t, y, p)
    dydt = zeros(4, 1);
    dydt(2) = (...
        + p.glucagon_blood_pancreas * y(2) * p.V_blood ...
        - p.glucagon_pancreas_portal * y(2) * p.V_pancreas ...
    ) / p.V_pancreas;
end

function dydt = somatostatin(t, y, p)
    dydt = zeros(4, 1); 
    dydt(3) = (...
        + p.somatostatin_blood_pancreas * y(3) * p.V_blood ...
        - p.somatostatin_pancreas_portal * y(3) * p.V_pancreas ...
    ) / p.V_pancreas;
end

function dydt = glucose(t, y, p)
    dydt = zeros(4, 1);
    dydt(4) = (...
        + p.glucose_blood_pancreas * y(4) * p.V_blood ...
        - p.glucose_pancreas_portal * y(4) * p.V_pancreas ...
    ) / p.V_pancreas;
end

