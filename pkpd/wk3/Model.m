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
            obj.params = struct();
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
        
        function r2 = R2(sample_data, fitted_data)
            arguments (Input)
                sample_data (1, :) double
                fitted_data (1, :) double
            end
            arguments (Output)
                r2 double
            end
            SS_res = sum((sample_data - fitted_data).^2);
            SS_tot = sum((sample_data - mean(sample_data)).^2);
            r2 = 1 - (SS_res/SS_tot);
        end

        function [t, y] = euler(odefun, tspan, y0, h)
            arguments
                odefun
                tspan (1, 2)
                y0 (1, :)
                h double = 0.01
            end
            % odefun: Function handle for the ODE (dy/dt = f(t, y))
            % tspan: Time span [t0 tf]
            % y0: Initial condition
            % h: Time step size
            
            % Number of time steps
            numSteps = ceil((tspan(2) - tspan(1)) / h);
            
            % Preallocate arrays for time and solution
            t = zeros(numSteps + 1, 1);
            y = zeros(numSteps + 1, length(y0));
            
            % Set initial conditions
            t(1) = tspan(1);
            y(1, :) = y0;
            
            % Euler integration loop
            for i = 1:numSteps
                % Compute the next value using Euler's method
                dydt = odefun(t(i), y(i, :));
                y(i + 1, :) = y(i, :)' + h .* dydt;
                t(i + 1) = t(i) + h;
            end
        end

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
            dydt(2) = p.ka*y(1) - p.kcl*y(2);
            dydt(3) = p.kcl*y(2);
        end
    end
end