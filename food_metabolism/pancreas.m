% TODO : minor optimizations involving passing by reference to decrease
% overhead from dydt initialization and hopefully a static inline to reduce
% function overhead (need to figure out how to compile this a c code)

% TODO : insulin needs to be self regulating
% TODO : glucagon needs to be self regulating
% TODO : somatostatin needs to be self regulating
% TODO : redefine glucose to modulate insulin, glucagon, and somatostatin secretion

classdef Pancreas
    methods (Static)

        function dydt = pancreas(...
            t, ...
            Ins, RIns, ...
            Glc, RGlc,...
            Somat, RSomat, ...
            Glu, RGlu, ...
            p_pancreas, ...
        )
            arguments (Input)
                t          double % timestamp
                Ins        double % Free Insulin Blood Conc.
                Rins       double % Available Insulin Receptor Conc.
                Glc        double % Free Glucagon Blood Conc. 
                RGlc       double % Available Glucagon Receptor Conc.
                Somat      double % Free Somatostatin Blood Conc.
                Rsomat     double % Available Somatostatin Receptor Conc.
                Glu        double % Free Glucose Blood Conc.
                RGlu       double % Available Glucose Receptor Conc.
                p_pancreas struct % parameter struct
            end
            arguments (Output)
                dydt (1, :) double
            end
            dinsulin = insulin(Ins, RIns, p_pancreas);
            dglucagon = glucagon(Glc, RGlc, p_pancreas);
            dsomatostatin = somatostatin(Somat, RSomat, p_pancreas);
            dglucose = glucose(Glu, RGlu, p_pancreas);

            dydt = [
                dinsulin, ...
                dglucagon, ...
                dsomatostatin, ...
                dglucose, ...
            ];
        end

        function dinsulin = insulin(I, R, IR, p)
            arguments (Input)
                % t double % timestamp
                I  double % Free Blood Insulin concentration
                R  double % Free Insulin Receptor concentration
                IR double % Bound Insulin
                p  struct % Insulin parameter struct
            end
            arguments (Output)
                dinsulin (1, 3) double
            end
            [dI, dR, dIR] = Ligand.ligandODE(I, R, IR, p.kinsulin_binding_on, p.kinsulin_binding_off);
            dinsulin = dinsulin + Ligand.inhibitory_secretion(...
                I, R, IR, ...
                p.kinsulin_binding_on, p.kinsulin_binding_off, ...
                p.kinsulin_baseline_secretion_inhibitory, ...
                p.kinsulin_tuner, ...
             );
            % dinsulin = dinsulin + glucose % upregulate
            % dinsulin = dinsulin + glucagon % inhibit
            % dinsulin = dinsulin + somatostatin % inhibit
            % dinsulin = dinsulin + proteins % upregulate
            % dinsulin = dinsulin + fatty acids % upregulate ? 
        end


        function dydt = glucagon(t, G, R, p)
            arguments (Input)
                t double
                G double
                R double
                p struct
            end
            arguments (Output)
                dglucose double
            end
            % L = Ligand.ligand_bound(I, R, p.kglucagon_on, pkglucagon_off);

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
