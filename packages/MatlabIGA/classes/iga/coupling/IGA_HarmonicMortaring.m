classdef IGA_HarmonicMortaring < DAE_Element
    properties (Access = public)
        MaterialParameter = 1/(4*pi*1e-7);
        Region1, Region2;
        CouplingInterface1, CouplingInterface2;
        a1DOFs;
        a2DOFs;
        couplingDOFs;
        RotationFunction;
    end

    methods (Access = public)
        function obj = IGA_HarmonicMortaring(region1, sides1, region2, sides2, materialParameter, nGaussPoints)
            if exist("materialParameter", "var") && ~isempty(materialParameter)
                obj.MaterialParameter = materialParameter;
            end
            if ~exist("nGaussPoints", "var")
                obj.CouplingInterface1 = IGA_CouplingInterface(region1, sides1);
                obj.CouplingInterface2 = IGA_CouplingInterface(region2, sides2);
            else
                obj.CouplingInterface1 = IGA_CouplingInterface(region1, sides1, nGaussPoints);
                obj.CouplingInterface2 = IGA_CouplingInterface(region2, sides2, nGaussPoints);
            end
            assert(abs(obj.CouplingInterface1.Radius - obj.CouplingInterface2.Radius) < 1e-10, "Coupling interfaces must have the same radius for Mortaring!")
            
            obj.Region1 = region1;
            obj.Region2 = region2;
            obj.addExternalElement(region1)
            obj.addExternalElement(region2)
        end

        function setHarmonics(obj, sinValues, cosValues)
            obj.CouplingInterface1.setHarmonics(sinValues, cosValues);
            obj.CouplingInterface2.setHarmonics(sinValues, cosValues);
            obj.NumberDOF = numel(obj.CouplingInterface1.HarmonicsSin) + numel(obj.CouplingInterface1.HarmonicsCos);
            obj.a1DOFs = 1:obj.Region1.NumberDOF;
            obj.a2DOFs = obj.Region1.NumberDOF+(1:obj.Region2.NumberDOF);
            obj.couplingDOFs = obj.NumberDOFext+(1:obj.NumberDOF);
        end

        function setRotationAngle(obj, angleRad)
            obj.CouplingInterface2.setRotationAngle(angleRad);
            % rotate region
            obj.Region1.RotationAngle = angleRad;
        end

        function [Mnonlin, Mlin] = MassMatrix(obj, t, y)
            Mnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            if nargout == 2
                Mlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            end
        end

        function [Knonlin, Klin] = StiffnessMatrix(obj, t, y)
            if ~isempty(obj.RotationFunction)
                obj.setRotationAngle(obj.RotationFunction(t))
            end
            % Nonlinear because matrix can change with rotation angle
            Knonlin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
            Knonlin(obj.a2DOFs, obj.couplingDOFs) = -obj.CouplingInterface2.CouplingMatrix * obj.MaterialParameter;
            Knonlin(obj.couplingDOFs, obj.a2DOFs) = -obj.CouplingInterface2.CouplingMatrix' * obj.MaterialParameter;
            if nargout == 2
                Klin = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
                Klin(obj.a1DOFs, obj.couplingDOFs) = obj.CouplingInterface1.CouplingMatrix * obj.MaterialParameter;
                Klin(obj.couplingDOFs, obj.a1DOFs) = obj.CouplingInterface1.CouplingMatrix' * obj.MaterialParameter;
            end
        end

        function [Fnonlin, Flin] = ForceVector(obj, t, y)
            Fnonlin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            if nargout == 2
                Flin = sparse(obj.NumberDOFext+obj.NumberDOF, 1);
            end
        end

        function J = JacobianMatrix(obj, t, y)
            J = sparse(obj.NumberDOFext+obj.NumberDOF, obj.NumberDOFext+obj.NumberDOF);
        end

        function postprocess(obj, ts, ys)
            a1 = ys(obj.a1DOFs, :);
            a2 = ys(obj.a2DOFs, :);
            obj.Solution.Time = ts;
            obj.Solution.Values = ys(obj.couplingDOFs, :);
            for iT = 1:numel(ts)
                t = ts(iT);
                if ~isempty(obj.RotationFunction)
                    obj.setRotationAngle(obj.RotationFunction(t));
                end
                obj.Solution.Torque(iT) = -obj.MaterialParameter*a2(:,iT)'*obj.CouplingInterface2.CouplingMatrixDer*obj.Solution.Values(:,iT);
            end
        end

        function T = calcTorqueBrBtRt(obj)
            T = 0;
            a = obj.Region1.MagneticPotential;
            for iBnd = 1:numel(obj.CouplingInterface1.Boundaries)
                patchDOFs = obj.Region1.Spaces.gnum{obj.Region1.Boundaries(obj.CouplingInterface1.Boundaries(iBnd)).patches};
                A_z = a(patchDOFs)/obj.Region1.Length;
                T = T + obj.MaterialParameter * obj.Region1.Length * ...
                op_Br_Bt_dGamma(obj.CouplingInterface1.SpacesEval{iBnd}, obj.CouplingInterface1.MeshesEval{iBnd}, A_z);
            end
        end

        function T = calcTorqueBrBtSt(obj)
            T = 0;
            a = obj.Region2.MagneticPotential;
            for iBnd = 1:numel(obj.CouplingInterface2.Boundaries)
                patchDOFs = obj.Region2.Spaces.gnum{obj.Region2.Boundaries(obj.CouplingInterface2.Boundaries(iBnd)).patches};
                A_z = a(patchDOFs)/obj.Region2.Length;
                T = T + obj.MaterialParameter * obj.Region2.Length * ...
                op_Br_Bt_dGamma(obj.CouplingInterface2.SpacesEval{iBnd}, obj.CouplingInterface2.MeshesEval{iBnd}, A_z);
            end
        end

        function [q1, q2] = calcHeatFlux(obj)
            q1 = 0;
            T = obj.Region1.Temperature;
            for iBnd = 1:numel(obj.CouplingInterface1.Boundaries)
                patchDOFs = obj.Region1.Spaces.gnum{obj.Region1.Boundaries(obj.CouplingInterface1.Boundaries(iBnd)).patches};
                t_loc = T(patchDOFs);
                q1 = q1 + obj.MaterialParameter * ...
                op_gradTn_dGamma(obj.CouplingInterface1.SpacesEval{iBnd}, obj.CouplingInterface1.MeshesEval{iBnd}, t_loc);
            end
            q2 = 0;
            T = obj.Region2.Temperature;
            for iBnd = 1:numel(obj.CouplingInterface2.Boundaries2)
                patchDOFs = obj.Region2.Spaces.gnum{obj.Region2.Boundaries(obj.CouplingInterface2.Boundaries(iBnd)).patches};
                t_loc = T(patchDOFs);
                q2 = q2 + obj.MaterialParameter * ...
                op_gradTn_dGamma(obj.CouplingInterface2.SpacesEval{iBnd}, obj.CouplingInterface2.MeshesEval{iBnd}, t_loc);
            end
        end

        function q = calcHeatFluxMortar(obj)
            q = sum(obj.CouplingInterface1.CouplingMatrix, 1)*obj.Solution.Values;
        end

        function plotCouplingValues(obj)
            figure;
            idxsin = mod(obj.CouplingInterface1.CouplingIndices, 2) == 1;
            idxcos = mod(obj.CouplingInterface1.CouplingIndices, 2) == 0;
            sinValues = obj.Solution.Values(idxsin);
            cosValues = obj.Solution.Values(idxcos);
            xlims = [min(obj.CouplingInterface1.HarmonicsAll)-1, max(obj.CouplingInterface1.HarmonicsAll)+1];
            subplot(2,1,1)
            bar(obj.CouplingInterface1.HarmonicsSin, sinValues);
            xticks(obj.CouplingInterface1.HarmonicsSin);
            xlim(xlims);
            title("Sin Values");
            subplot(2,1,2)
            bar(obj.CouplingInterface1.HarmonicsCos, cosValues);
            xticks(obj.CouplingInterface1.HarmonicsCos);
            xlim(xlims);
            title("Cos Values");
        end
    end
end