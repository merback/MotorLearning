function dfdP = op_D_DP_gradu_nu_gradv_v_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, v, dCdP, mat, len, patch_list)

    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    nParams = size(dCdP, 3);
    
    for iPatch = patch_list
        [r, c, vals] = op_D_DP_gradu_nu_gradv_v(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), v(spv.gnum{iPatch}), dCdP(spg.gnum{iPatch},:,:), mat, len);
    
        r = spu.gnum{iPatch}(r);
    
        if (~isempty (spv.dofs_ornt))
            vals = spv.dofs_ornt{iPatch}(r)' .* vals;
        end
    
        rows{iPatch} = r;
        cols{iPatch} = c;
        values{iPatch} = vals;
    end
    
    dfdP = sparse(cell2mat(rows), cell2mat(cols), cell2mat(values), spu.ndof, nParams);
end