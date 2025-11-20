classdef MNA_Inductance < DAE_Element
    properties
        L;
    end
    methods
        function obj = MNA_Inductance(L)
            obj.NumberDOF = 2;
            obj.L = L;
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % d/dt Phi = AL'*phi = uL
        % 0 = Phi - Li;
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(4,4);
            Mlin = [0 0 0 0; 
	                0 0 0 0; 
	                0 0 1 0; 
	                0 0 0 0];
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(4,4);
            Klin = [0 0 0 1; 
                    0 0 0 -1; 
                    1 -1 0 0; 
                    0 0 1 -obj.L];
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(4,1);
            Flin = sparse(4,1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(4,4);
        end

        function postprocess(obj, ts, ys)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys(1:2, :);
            obj.Solution.Flux = ys(3, :);
            obj.Solution.Current = ys(4, :);
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
            subplot(3,1,2)
            plot(obj.Solution.Time, obj.Solution.Flux);
            ylabel("Flux (Wb)")
            xlabel("Time (s)")
            title("Flux")
            grid on
            subplot(3,1,3)
            plot(obj.Solution.Time, obj.Solution.Current);
            ylabel("Current (A)")
            xlabel("Time (s)")
            title("Current")
            grid on
        end
    end
end

