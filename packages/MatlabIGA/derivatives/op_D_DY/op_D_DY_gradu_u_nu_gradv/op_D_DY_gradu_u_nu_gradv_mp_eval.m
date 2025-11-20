function A = op_D_DY_gradu_u_nu_gradv_mp_eval(spu, spuEval, spv, spvEval, msh, mshEval, u, material, patch_list)

    rows = cell(msh.npatch, 1);
    cols = cell(msh.npatch, 1);
    values = cell(msh.npatch, 1);
    
    for iptc = patch_list
        [r, c, vals] = op_D_DY_gradu_u_nu_gradv(spuEval{iptc}, spvEval{iptc}, mshEval{iptc}, u(spu.gnum{iptc}), material);
        rows{iptc} = spv.gnum{iptc}(r);
        cols{iptc} = spu.gnum{iptc}(c);
    
        if (~isempty (spv.dofs_ornt))
            vals = spv.dofs_ornt{iptc}(r)' .* vals;
        end
        if (~isempty (spu.dofs_ornt))
            vals = vals .* spu.dofs_ornt{iptc}(c)';
        end
        values{iptc} = vals;
    end
    
    A = sparse(cell2mat(rows), cell2mat(cols), cell2mat(values), spv.ndof, spu.ndof);
end