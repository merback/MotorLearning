classdef IGA_EXCmagnet < IGA_EXC
    properties (Access = public)
        Material;
    end

    methods
        function obj = IGA_EXCmagnet(parentRegion, material, patches)
            obj@IGA_EXC(parentRegion);
            obj.IsConstant = true;
            obj.Material = material;
            obj.Patches = patches;
            obj.updateExcitation();
        end

        function f = getEXCvalues(obj, t)
            f = obj.RHS;
        end

        function updateExcitation(obj)
            obj.RHS = op_gradv_Br_bot_mp(obj.ParentRegion.Spaces, obj.ParentRegion.Meshes, obj.Material, obj.Patches);
        end
    end
end