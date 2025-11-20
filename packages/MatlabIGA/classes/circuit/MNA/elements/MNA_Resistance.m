classdef MNA_Resistance < DAE_Element
    properties
        Conductivity;
        Length;
        CrossSection;

        Resistance;
        Conductance;
    end
    methods
        function obj = MNA_Resistance(varargin)
            obj.NumberDOF = 0;
            if nargin == 1
                obj.Resistance = varargin{1};
            elseif nargin == 3
                obj.Conductivity = varargin{1};
                obj.Length = varargin{2};
                obj.CrossSection = varargin{3};
                obj.Resistance = obj.CrossSection/(obj.Conductivity*obj.Length);
            else
                error([class(obj) ': Wrong number of input arguments!'])
            end
            obj.Conductance = 1/obj.Resistance;
        end

        %%%%%%%% element equations
        % 0 = i = AR*G*AR'*phi
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(2, 2);
            Mlin = sparse(2, 2);
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(2, 2);
            if nargout == 2
                Klin = obj.Conductance* [1, -1; -1, 1];
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(2,1);
            Flin = sparse(2,1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(2,2);
        end

        %%%%%%%% postprocessing
        function postprocess(obj, ts, ys)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = [1, -1]*ys;
            obj.Solution.Current = obj.Conductance*obj.Solution.Voltage;
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

