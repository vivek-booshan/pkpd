% TODO: Need to verify if classdef overhead is worth wrapping in static class
% 
classdef Ligand
    methods (Static)

        function LR = ligand_bound(L, R, kon, koff)
            % Assumes Pseudo steady state assumption
            arguments (Input)
                L    double % Free Ligand Concentration
                R    double % Free Receptor Concentration
                kon  double % L + R --> LR
                koff double % LR --> L + R
            end
            arguments (Output)
                LR double % Bound Ligand Concentration
            end

            LR = kon * L * R / koff
        end

        function Kd = ligand_Kd(L, R, LR)
            arguments (Input)
                L double % Free Ligand Concentration
                R double % Free Receptor Concentration
                LR double % Bound Ligand Concentration
            end
            arguments (Output)
                Kd double % Binding Dissociation constant
            end
            Kd = L * R / LR
        end

        function f = fractional_occupancy_Kd(L, Kd)
            f = L / (L + Kd)
        end

        function f = fractional_occupancy_LR(R, LR)
            f = 1 / (1 + (R / LR))
        end

        function f = fractional_inhibitory_effect(f_occupancy, k)
            denom = 1 + exp(k * (f - 0.5))
            f = 1/denom
        end

        % TODO : what even is baseline secretion? Need to redefine in terms of max secretion
        function keff = inhibitory_secretion(L, R, LR, kon, koff, kbaseline_secretion, k)
            Kd = L * R / LR;
            fo = L / (L + Kd);
            f = 1 / (1 + exp(k * (f - 0.5)));
            keff = kbaseline_secretion * f;
        end

        function dligand = ligandODE(L, R, LR, kon, koff)
            dL  = -kon * L * R + koff * LR;
            dR  = -kon * L * R + koff * LR;
            dLR = +kon * L * R - koff * LR;
        end

end
