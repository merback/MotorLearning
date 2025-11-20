classdef MNA_InductanceMutual2 < DAE_Element
    properties
        ts;
        Ls;
        AL;
    end
    methods
        function obj = MNA_InductanceMutual2()
            obj.NumberDOF = 2;
            obj.AL = sparse([1:2:2*obj.NumberDOF, 2:2:2*obj.NumberDOF], [1:obj.NumberDOF, 1:obj.NumberDOF], ...
                    [ones(1, obj.NumberDOF), -ones(1, obj.NumberDOF)], 2*obj.NumberDOF, obj.NumberDOF);
        end

        function setL(obj, t, L)
            L = reshape(L, 1, []);
            if any(obj.ts == t)
                ind = find(obj.ts==t);
                obj.Ls(ind,:) = L;
            else
                obj.ts = [obj.ts, t];
                obj.Ls = [obj.Ls; L];
                [obj.ts, order] = sort(obj.ts);
                obj.Ls = obj.Ls(order,:);
            end
        end

        %%%%%%%% element equations
        % 0 = AL*iL
        % L d/dt i = AL'*phi
        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
            L = interp1(obj.ts, full(obj.Ls), t, "linear", "extrap");
            L = reshape(L, 2, 2);
            Mnonlin(2*obj.NumberDOF+(1:obj.NumberDOF),2*obj.NumberDOF+(1:obj.NumberDOF)) = L;
            if nargout == 2
                Mlin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            Knonlin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
            if nargout == 2
                Klin = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
                Klin(1:2*obj.NumberDOF, 2*obj.NumberDOF + (1:obj.NumberDOF)) = obj.AL; 
                Klin(2*obj.NumberDOF + (1:obj.NumberDOF), 1:2*obj.NumberDOF) = obj.AL'; 
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(3*obj.NumberDOF, 1);
            Flin = sparse(3*obj.NumberDOF, 1);
        end
        
        function J = JacobianMatrix(obj, t, y)
            J = sparse(3*obj.NumberDOF, 3*obj.NumberDOF);
        end

        function postprocess(obj, ts, ys)
            obj.Solution.Time = ts;
            obj.Solution.Voltage = obj.AL'*ys(1:2*obj.NumberDOF, :);
            obj.Solution.Current = ys(2*obj.NumberDOF+(1:obj.NumberDOF), :);
            for ti = 1:numel(ts)
                L = interp1(obj.ts, full(obj.Ls), ts(ti), "linear", "extrap");
                L = reshape(L, 2, 2);
                obj.Solution.Flux(:,ti) = L*obj.Solution.Current(:,ti);
            end
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

