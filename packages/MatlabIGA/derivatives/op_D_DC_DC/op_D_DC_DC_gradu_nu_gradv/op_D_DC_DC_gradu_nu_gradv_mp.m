
function dKdCdC = op_D_DC_DC_gradu_nu_gradv_mp(spu, spv, spg, msh, u, material, patch_list)
    
    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list

         [ind, vals] = op_D_DC_DC_gradu_nu_gradv_tp (spu.sp_patch{iPatch}, spv.sp_patch{iPatch}, spg.sp_patch{iPatch}, msh.msh_patch{iPatch}, u(spu.gnum{iPatch}), material);

        ind(:,1) = spv.gnum{iPatch}(ind(:,1));
        ind(:,2) = spu.gnum{iPatch}(ind(:,2));
        ind(:,3) = spg.gnum{iPatch}(ind(:,3));
        ind(:,5) = spg.gnum{iPatch}(ind(:,5));
    
        if (~isempty (spu.dofs_ornt))
            vs = spu.dofs_ornt{iPatch}(ind(:,2))' .* vs;
        end
        if (~isempty (spv.dofs_ornt))
            vs = vs .* spv.dofs_ornt{iPatch}(ind(:,1))';
        end
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DC_DC_gradu_nu_gradv_mp: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iPatch}(ind(:,3))';
        end
    
        indices{iPatch} = ind;
        values{iPatch} = vals;
    end

    dKdCdC = sptensor (cell2mat(indices), cell2mat(values), [spv.ndof, spu.ndof, spg.ndof, msh.ndim, spg.ndof, msh.ndim]);
end
