% TODO : minor optimizations involving passing by reference to decrease
% overhead from dydt initialization and hopefully a static inline to reduce
% function overhead (need to figure out how to compile this a c code)

% TODO : insulin needs to be self regulating
% TODO : glucagon needs to be self regulating
% TODO : somatostatin needs to be self regulating
% TODO : redefine glucose to modulate insulin, glucagon, and somatostatin secretion

classdef Pancreas
    methods (Static)

        function dydt = pancreas(t, y, p)
            dydt = zeros(4, 1);
            dinsulin = insulin(t, y, p);
            dglucagon = glucagon(t, y, p);
            dsomatostatin = somatostatin(t, y, p);
            dglucose = glucose(t, y, p);

            % dydt = dinsulin + dglucagon + dsomatostatin + dglucose;
        end

        function dinsulin = insulin(t, I, R, p)
            arguments (Input)
                t double % timestamp
                I double % Free Blood Insulin concentration
                R double % Free Insulin Receptor concentration
                p struct % Insulin parameter struct
            end
            arguments (Output)
                dinsulin double
            end

            Kd = Ligand.ligand_Kd(L, R, LR);
            fo = Ligand.fractional_occupancy(L, Kd);
            dinsulin = p.k_insulin_secretion_insulin * ligand.fractional_inhibitory_effect(fo, 15) * I;
            % dinsulin = dinsulin + glucose % upregulate
            % dinsulin = dinsulin + glucagon % inhibit
            % dinsulin = dinsulin + somatostatin % inhibit
            % dinsulin = dinsulin + proteins % upregulate
            % dinsulin = dinsulin + fatty acids % upregulate ? 
        end


        function dydt = glucagon(t, y, p)
            dglucagon = (...
                + p.glucagon_blood_pancreas * y(2) * p.V_blood ...
                - p.glucagon_pancreas_portal * y(2) * p.V_pancreas ...
            ) / p.V_pancreas;
        end

        function dydt = somatostatin(t, y, p)
            dsomatostatin = (...
                + p.somatostatin_blood_pancreas * y(3) * p.V_blood ...
                - p.somatostatin_pancreas_portal * y(3) * p.V_pancreas ...
            ) / p.V_pancreas;
        end

        function dydt = glucose(t, y, p)
            dglucose = (...
                + p.glucose_blood_pancreas * y(4) * p.V_blood ...
                - p.glucose_pancreas_portal * y(4) * p.V_pancreas ...
            ) / p.V_pancreas;
        end

    end
end
