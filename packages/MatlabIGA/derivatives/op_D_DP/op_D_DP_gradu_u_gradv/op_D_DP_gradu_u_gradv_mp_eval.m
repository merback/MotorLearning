function dKudP = op_D_DP_gradu_u_gradv_mp_eval(spu, spuEval, spv, spvEval, spg, spgEval, msh, mshEval, u, dCdP, patch_list)

    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    vals = cell(msh.npatch, 1);
    nParams = size(dCdP, 3);
    
    for iPatch = patch_list
        [r, c, v] = op_D_DP_gradu_u_gradv(spuEval{iPatch}, spvEval{iPatch}, spgEval{iPatch}, mshEval{iPatch}, u(spu.gnum{iPatch}), dCdP(spg.gnum{iPatch},:,:));
    
        r = spv.gnum{iPatch}(r);
    
        if (~isempty (spv.dofs_ornt))
            v = spv.dofs_ornt{iPatch}(r(:,1))' .* v;
            warning("DOF ornts for u not yet checked");
        end
    
        rows{iPatch} = r;
        cols{iPatch} = c;
        vals{iPatch} = v;
    end
    
    dKudP = sparse(cell2mat(rows), cell2mat(cols), cell2mat(vals), spv.ndof, nParams);
end