classdef BC_AntiPeriodic < BC
    properties (Access = public)
        LeftDOFs, 
        RightDOFs;
    end

    methods (Access = public)
        function obj = BC_AntiPeriodic(parentElement, leftDOFs, rightDOFs)
            leftDOFs = reshape(leftDOFs, 1, []); rightDOFs = reshape(rightDOFs, 1, []);
            assert(isempty(intersect(leftDOFs, rightDOFs)), "Left and right DOFs must not be identical!");
            % Reduce DOFs for possible multiplicity with other Boundary Conditions
            obj@BC(parentElement, union(leftDOFs, rightDOFs, 'stable'));
            obj.LeftDOFs = intersect(leftDOFs, obj.BoundaryDOFs, 'stable');
            obj.RightDOFs = intersect(rightDOFs, obj.BoundaryDOFs, 'stable');
            % only set LeftDOFs as the DOFs which should be removed
            obj.RemoveDOFs = obj.LeftDOFs;
            % all others are independent
            obj.KeepDOFs = setdiff((1:obj.ParentElement.NumberDOF), obj.RemoveDOFs);
            obj.Name = "AntiPeriodic";

            C = sparse(repmat(1:numel(obj.LeftDOFs), 1, 2), [obj.LeftDOFs, obj.RightDOFs], ...
                ones(1, numel(obj.LeftDOFs)+numel(obj.RightDOFs)), ...
                numel(obj.RemoveDOFs), numel(obj.RemoveDOFs)+numel(obj.KeepDOFs));
            obj.CK = C(:, obj.KeepDOFs);
            obj.CR = C(:, obj.RemoveDOFs);
        end

        function [b, db]  = getBCvalues(obj, t)
            b = sparse(numel(obj.RemoveDOFs), 1);
            db = sparse(numel(obj.RemoveDOFs), 1);
        end
    end
end