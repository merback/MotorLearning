function varargout = op_D_DY_gradu_u_nu_gradv_tp(space1, space2, msh, u, material)

    for idim = 1:msh.ndim
        size1 = size (space1.sp_univ(idim).connectivity);
        size2 = size (space2.sp_univ(idim).connectivity);
        if ((size1(2) ~= size2(2)) || (size2(2) ~= msh.nel_dir(idim)))
            error ('One of the discrete spaces is not associated to the mesh')
        end
    end
    
    A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);
    
    for iel = 1:msh.nel_dir(1)
        msh_col = msh_evaluate_col (msh, iel);
        sp1_col = sp_evaluate_col(space1, msh_col, 'value', false, 'gradient', true);
        sp2_col = sp_evaluate_col(space2, msh_col, 'value', false, 'gradient', true);
    
        A = A + op_D_DY_gradu_u_nu_gradv(sp1_col, sp2_col, msh_col, u, material);
    end
    
    if (nargout == 1)
        varargout{1} = A;
    elseif (nargout == 3)
        [rws, cls, vals] = find (A);
        varargout{1} = rws;
        varargout{2} = cls;
        varargout{3} = vals;
    end
end