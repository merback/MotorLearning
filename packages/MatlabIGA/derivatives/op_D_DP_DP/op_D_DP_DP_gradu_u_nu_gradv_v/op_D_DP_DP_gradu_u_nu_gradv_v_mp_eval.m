
function dKdPdP = op_D_DP_DP_gradu_u_nu_gradv_v_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, v, dCdP, d2CdP2, mat, len, patch_list)
    
    nParams = size(dCdP, 3);

    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list

        [r, c, vals] = op_D_DP_DP_gradu_u_nu_gradv_v(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), v(spv.gnum{iPatch}), dCdP(spg.gnum{iPatch},:,:), d2CdP2(spg.gnum{iPatch},:,:,:), mat, len);
    
        rows{iPatch} = r;
        cols{iPatch} = c;
        values{iPatch} = vals;
    end

    dKdPdP = sparse (cell2mat(rows), cell2mat(cols), cell2mat(values), nParams, nParams);
end
