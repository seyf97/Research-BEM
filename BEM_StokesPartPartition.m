function assembly = BEM_StokesPartPartition(mesh, assembly)
    % Partitioning the system matrices to to form the linear system
    % required for easily reaching the solution
    
    
    %%
    %{
    ***********************************************************************
        Reforming System Matrices
    ***********************************************************************
    %}
    if ~mesh.findForce
        %
        assembly.A = assembly.HPP - mtimes(assembly.GP0, mtimes(assembly.G00Inv, assembly.H0P));
    end
    %
    assembly.B = assembly.GPP - mtimes(assembly.GP0, mtimes(assembly.G00Inv, assembly.G0P));
    
    %
    assembly.C = assembly.HP0 - mtimes(assembly.GP0, mtimes(assembly.G00Inv, assembly.H00));
    
    
    %%
    %{
    ***********************************************************************
        Inversion of B Matrix
    ***********************************************************************
    %}
    
    assembly.BInv = eye(2 * mesh.N_E_P_O) / assembly.B;
    
    
    %%
    %{
    ***********************************************************************
        Finalizing Repartition
    ***********************************************************************
    %}
    if ~mesh.findForce
        %
        assembly.T = mtimes(assembly.F, mtimes(assembly.BInv, mtimes(assembly.A, assembly.M)));
    end
    
    %
    assembly.R = mtimes(assembly.F, mtimes(assembly.BInv, assembly.C));
end