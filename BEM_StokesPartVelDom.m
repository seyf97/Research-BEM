function [assembly, values] = BEM_StokesPartVelDom(mesh, assembly, values)
    % Calculation of the velocity field inside the computational domain
    
    
    %%
    %{
    ***********************************************************************
        System Matrices 
    ***********************************************************************
    
    Matrices necessary for generating the velocity field in the
    domain
    %}
    assembly = BEM_StokesPartVelDomAssem(mesh, assembly);
    
    %%
    %{
    Velocity Field inside the Computational Domain
    %}
    values.uD = assembly.W * values.t - assembly.P * values.u;
end