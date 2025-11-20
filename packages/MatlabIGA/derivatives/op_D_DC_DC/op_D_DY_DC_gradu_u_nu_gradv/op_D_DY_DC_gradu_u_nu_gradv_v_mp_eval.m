
function dKdYdC = op_D_DY_DC_gradu_u_nu_gradv_v_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, v, material, len, patch_list)
    
    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list

         [ind, vals] = op_D_DY_DC_gradu_u_nu_gradv_v(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), v(spv.gnum{iPatch}), material, len);

        ind(:,1) = spu.gnum{iPatch}(ind(:,1));
        ind(:,2) = spg.gnum{iPatch}(ind(:,2));
    
        if (~isempty (spv.dofs_ornt))
            vs = vs .* spv.dofs_ornt{iPatch}(ind(:,1))';
        end
        if (~isempty (spg.dofs_ornt))
            warning("op_D_DY_DC_gradu_u_nu_gradv_v_mp_eval: dofs_ornt not tested yet")
            vs = vs .* spg.dofs_ornt{iPatch}(ind(:,2))';
        end
    
        indices{iPatch} = ind;
        values{iPatch} = vals;
    end

    dKdYdC = sptensor (cell2mat(indices), cell2mat(values), [spu.ndof, spg.ndof, msh.ndim]);
end
