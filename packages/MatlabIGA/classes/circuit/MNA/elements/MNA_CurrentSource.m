classdef MNA_CurrentSource < DAE_Element
    properties
        Current;
    end
    methods
        function obj = MNA_CurrentSource(ic)
            obj.NumberDOF = 0;
            if nargin == 1
                if isnumeric(ic)
                    obj.Current = @(t) ic*ones(size(t));
                elseif isa(ic, 'function_handle')
                    obj.Current = ic;
                else
                    error([class(obj) ': Wrong input type!'])
                end
            else
                error([class(obj) ': Wrong number of input arguments!'])
            end
        end

        %%%%%%%% element equations
        % 0 = AI*iI
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2,2);
            Mlin = sparse(2,2);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2,2);
            Klin = sparse(2,2);
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = [1; -1]*obj.Current(t);
            Flin = sparse(2,1);
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(2,2);
        end

        function postprocess(obj, ts, ys)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys(1:2,:);
            obj.Solution.Current = obj.Current(ts);
        end

        function plotSolution(obj)
            figure()
            sgt = sgtitle(obj.Name);
            sgt.Interpreter = "none";
            subplot(2,1,1)
            plot(obj.Solution.Time, obj.Solution.Voltage)
            ylabel("Voltage (V)")
            xlabel("Time (s)")
            title("Voltage")
            grid on
            subplot(2,1,2)
            plot(obj.Solution.Time, obj.Solution.Current);
            ylabel("Current (A)")
            xlabel("Time (s)")
            title("Current")
            grid on
        end
    end
end

