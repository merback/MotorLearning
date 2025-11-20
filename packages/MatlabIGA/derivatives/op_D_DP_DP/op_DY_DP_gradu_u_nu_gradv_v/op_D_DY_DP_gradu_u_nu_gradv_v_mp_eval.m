
function dKdYdP = op_D_DY_DP_gradu_u_nu_gradv_v_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, v, dCdP, mat, len, patch_list)

    nParams = size(dCdP, 3);
    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list
    
        [r, c, vals] = op_D_DY_DP_gradu_u_nu_gradv_v(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), v(spv.gnum{iPatch}), dCdP(spg.gnum{iPatch},:,:), mat, len);
    
        r = spu.gnum{iPatch}(r);
    
        if (~isempty (spv.dofs_ornt))
            vs = vs .* spv.dofs_ornt{iPatch}(r)';
        end
    
        rows{iPatch} = r;
        cols{iPatch} = c;
        values{iPatch} = vals;
    end
    
    dKdYdP = sparse(cell2mat(rows), cell2mat(cols), cell2mat(values), spu.ndof, nParams);
end
