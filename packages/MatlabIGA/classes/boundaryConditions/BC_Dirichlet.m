classdef BC_Dirichlet < BC
    properties (Access = public)
        Value
        ValueDerivative
    end

    methods (Access = public)
        function obj = BC_Dirichlet(parentElement, boundaryDOFs, value, dvalue)
            obj@BC(parentElement, boundaryDOFs);
            obj.RemoveDOFs = obj.BoundaryDOFs;
            obj.KeepDOFs = setdiff(1:parentElement.NumberDOF, obj.RemoveDOFs);

            obj.Name = "Dirichlet";
            if ~exist("value", "var")
                value = 0;
            end
            if ~exist("dvalue", "var")
                dvalue = 0;
            end
            obj.ValueDerivative = dvalue;
            obj.Value = value;

            obj.CK = sparse(numel(obj.RemoveDOFs), numel(obj.KeepDOFs));
            obj.CR = sparse(eye(numel(obj.RemoveDOFs)));
        end

        function [b, db] = getBCvalues(obj, t)
            if isa(obj.Value, "function_handle")
                b = ones(numel(obj.RemoveDOFs),1)*obj.Value(t);
            else
                b = ones(numel(obj.RemoveDOFs),1)*obj.Value;
            end

            if isa(obj.ValueDerivative, "function_handle")
                db = ones(numel(obj.RemoveDOFs),1)*obj.ValueDerivative(t);
            else
                db = ones(numel(obj.RemoveDOFs),1)*obj.ValueDerivative;
            end
        end
    end
end