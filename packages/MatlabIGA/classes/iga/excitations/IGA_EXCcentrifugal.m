classdef IGA_EXCcentrifugal < IGA_EXC
    properties (Access = public)
        RotationSpeed;
    end

    methods
        function obj = IGA_EXCcentrifugal(parentRegion, rotationSpeed)
            obj@IGA_EXC(parentRegion);
            obj.Patches = 1:parentRegion.NumberPatches;
            obj.RotationSpeed = rotationSpeed;
            if isa(rotationSpeed, "function_handle")
                obj.IsConstant = false;
            else
                obj.IsConstant = true;
            end

            obj.updateExcitation()
        end

        function updateExcitation(obj)
            obj.RHS = zeros(obj.ParentRegion.NumberDOF, 1);
            for iMat = 1:numel(obj.ParentRegion.Materials)
                func = @(x,y) obj.ParentRegion.Materials(iMat).Material.getRho .* cat(1, reshape (cos(atan2(y,x)).* sqrt(x.^2 + y.^2), [1, size(x)]), reshape (sin(atan2(y,x)).* sqrt(x.^2 + y.^2), [1, size(x)]));
                obj.RHS = obj.RHS + op_f_v_mp(obj.ParentRegion.Spaces, obj.ParentRegion.Meshes, func , obj.ParentRegion.Materials(iMat).Patches);
            end
        end
        
        function f = getEXCvalues(obj, t)
            if obj.IsConstant
                omega = 2.*pi.*obj.RotationSpeed./60;
                f = omega^2*obj.RHS;
            else
                omega = 2.*pi.*obj.RotationSpeed(t)./60;
                f = omega^2*obj.RHS;
            end
        end
    end
end