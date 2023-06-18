function [assembly, values] = BEM_StokesPartPresDom(mesh, assembly, values)
    % Pressure Field in the Domain
    % 
    % 
    % mu : Dynamic Viscosity of the fluid [Pa.s]
    % 
    % N_E : Number of elements at the boundaries [1 x 1]
    % 
    % N_D : Number of points in the domain [1 x 1]
    % 
    % X_D : Matrix containing the positions of the nodes in the domain 
    %       [N_D x 2]
    % 
    % XM : Matrix containing the positions of the nodes at the boundary 
    %       [N_E x 2]
    % 
    % E_L : Vector containing the lengths of the elements at the boundary
    %       [N_E x 1]
    % 
    % E_N : Matrix containing the components of the outward unit normal 
    %       vectors of the elements at the boundary [N_E x 2]
    % 
    % Tg : Matrix containing the relative positions of the elements at the
    %      boundary [N_E x 2]
    % 
    % u : Matrix containing the velocity components of the elements at the
    %     boundary [2 * N_E x 1]
    % 
    % t : Matrix containing the traction components of the elements at the
    %     boundary [2 * N_E x 1]
    
    
    %%
    %{
    ***********************************************************************
        System Matrices 
    ***********************************************************************
    
    Matrices necessary for generating the pressure field in the
    domain
    %}
    assembly = BEM_StokesPartPresDomAssem(mesh, assembly);
    
    
    %%
    %{
    ***********************************************************************
    Pressure Field in the Domain
    ***********************************************************************
    %}
    values.pD = assembly.Q * values.t - assembly.X * values.u;
end