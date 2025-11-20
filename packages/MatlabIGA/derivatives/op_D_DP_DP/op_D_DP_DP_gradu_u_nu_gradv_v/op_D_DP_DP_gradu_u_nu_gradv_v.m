function varargout = op_D_DP_DP_gradu_u_nu_gradv_v(spu, spv, spg, msh, u, z, dCdP, d2CdP2, mat, len)

    % First order derivative
    dKdC = op_D_DC_gradu_u_nu_gradv_v(spu, spv, spg, msh, u, z, mat, len);
    
    % Second order derivative
    dKdCdC = op_D_DC_DC_gradu_u_nu_gradv_v(spu, spv, spg, msh, u, z, mat, len);
    
    % Product rule
    dKdPdP1 = ttt(sptensor(dKdC), sptensor(d2CdP2(:,1:msh.ndim,:,:)), [1,2], [1,2]);

    dKdPdP2 = ttt(sptensor(ttt(dKdCdC, sptensor(dCdP(:,1:msh.ndim,:)), [3,4], [1,2])), ...
        sptensor(dCdP(:,1:msh.ndim,:)), [1,2], [1,2]);

    dKdPdP = double(dKdPdP1 + dKdPdP2);

    if (nargout == 1 || nargout == 0)
        varargout{1} = dKdPdP;
    elseif (nargout == 3)
        [rows, cols, values] = find(dKdPdP);
        varargout{1} = rows;
        varargout{2} = cols;
        varargout{3} = values;
    else
        error ('op_D_DP_DP_gradu_u_nu_gradv_v: wrong number of output arguments')
    end
end
