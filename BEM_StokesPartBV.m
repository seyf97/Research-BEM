function values = BEM_StokesPartBV(mesh, assembly, values)
    % 
    
    
    %%
    %{
    ***********************************************************************
        Boundary Condition Implementation
    ***********************************************************************

    Rearranging the entries of the matrices according to the 
    boundary conditions to get them ready for being used in a 
    linear system the  solution of which would be the unknown 
    boundary values we are looking for
    %}
    [H_Imp, G_Imp] = BEM_StokesPartBC_Imp...
        (mesh.N_E, mesh.BC, assembly.H, assembly.G);
    
    
    %%
    %{
    ***********************************************************************
        Boudary Value Calculation
    ***********************************************************************

    Solving the linear system to get the unknown boundary 
    values
                                    &
    Generating the values of the velocity (both x and y 
    components) and  the traction (both x and y components) for
    each element

    The velocities at the channel nodes will be fully known after mesh.BC 
    Implementation
    %}
    
    % Solving for unknown velocity and traction values
    [values.u, values.t] =...
        BEM_StokesPartSolver(mesh.N_E, mesh.BC, H_Imp, G_Imp);
    
    %{
    Velocities at the channel nodes

    u0 : (2 * N0) x 1
    %}
    values.u0 = values.u(1 : 2 * mesh.N_E_0);
end