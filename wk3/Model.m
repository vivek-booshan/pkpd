classdef Model < handle
    properties
        params struct
        method logical = 0; % 0 if infusion, 1 if oral
    end

    methods
        function obj = Model(q, V, kcl, ka)
            arguments
                q double
                V double 
                kcl double 
                ka double = [] 
            end
            obj.params.q = q;
            obj.params.V = V;
            obj.params.kcl = kcl;

            if ~isempty(ka)
                obj.params.ka = ka;
                obj.method = 1;
            end
        end


    end
    methods (Static)
        
        function dydt = infusionODE(t, y, p)
            arguments
                t (1, :) double
                y (1, 2) double 
                p struct
            end
            dydt = zeros(2, 1);
            dydt(1) = p.q / p.V - p.kcl * y(1);
            dydt(2) = p.kcl * y(1);
        end

        function dydt = oralODE(t, y, p)
            arguments
                t (1, :) double
                y (1, 3) double
                p struct
            end
            dydt = zeros(3, 1);
            dydt(1) = -p.ka*y(1);
            dydt(2) = p.ka*y(1)/p.V - p.kcl*y(2);
            dydt(3) = p.kcl*y(2);
        end
    end
end