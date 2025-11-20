classdef MNA_Capacitance < DAE_Element
    properties
        C;
    end
    methods
        function obj = MNA_Capacitance(C)
            obj.NumberDOF = 0;
            obj.C = C;
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % d/dt Phi = AL'*phi = uL
        % 0 = Phi - Li;
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2,2);
            Mlin = [-obj.C ,obj.C; 
                    obj.C, -obj.C];
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2,2);
            Klin = sparse(2,2);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(2,1);
            Flin = sparse(2,1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(2,2);
        end

        function postprocess(obj, ts, ys)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys(1:2, :);
            % obj.Solution.Charge = ys(3, :); % TBD
            % obj.Solution.Current = ys(4, :); % TBD
        end

        function plotSolution(obj)
            figure()
            sgt = sgtitle(obj.Name);
            sgt.Interpreter = "none";
            subplot(3,1,1)
            plot(obj.Solution.Time, obj.Solution.Voltage)
            ylabel("Voltage (V)")
            xlabel("Time (s)")
            title("Voltage")
            grid on
            % subplot(3,1,2)
            % plot(obj.Solution.Time, obj.Solution.Flux);
            % ylabel("Flux (Wb)")
            % xlabel("Time (s)")
            % title("Flux")
            % grid on
            % subplot(3,1,3)
            % plot(obj.Solution.Time, obj.Solution.Current);
            % ylabel("Current (A)")
            % xlabel("Time (s)")
            % title("Current")
            % grid on
        end
    end
end

