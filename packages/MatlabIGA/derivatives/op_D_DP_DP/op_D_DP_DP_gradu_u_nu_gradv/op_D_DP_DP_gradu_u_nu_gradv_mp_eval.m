
function dKdPdP = op_D_DP_DP_gradu_u_nu_gradv_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, dCdP, d2CdP2, mat, patch_list)
    
    nParams = size(dCdP, 3);

    indices = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iPatch = patch_list

        [ind, vals] = op_D_DP_DP_gradu_u_nu_gradv(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), dCdP(spg.gnum{iPatch},:,:), d2CdP2(spg.gnum{iPatch},:,:,:), mat);

        ind(:,1) = spv.gnum{iPatch}(ind(:,1));
    
        indices{iPatch} = ind;
        values{iPatch} = vals;
    end

    dKdPdP = sptensor (cell2mat(indices), cell2mat(values), [spv.ndof, nParams, nParams]);
end
