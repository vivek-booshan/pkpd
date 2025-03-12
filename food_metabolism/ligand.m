% TODO: Need to verify if classdef overhead is worth wrapping in static class
% 
classdef Ligand
    methods (Static)

        function LR = ligand_bound(L, R, kon, koff)
            arguments
                L double % Free Ligand Concentration
                R double % Free Receptor Concentration
                kon double % L + R --> LR
                koff double % LR --> L + R
            end

            LR = kon * L * R / koff
        end

        function Kd = ligand_Kd(L, R, LR)
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

    end
end
