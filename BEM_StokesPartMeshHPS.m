function mesh = BEM_StokesPartMeshHPS(mesh)
    % Discretization of the channel wall
    % 
    %                         
    %                             Figure 1
    %   
    %
    %        <------------------------- Lx ------------------------->  
    %        ___________________                   __________________  
    %    ^  |                   \                 /                  | 
    %    |  |                    \               /                   | 
    %    |  |                     \             /                    | 
    %    |  |                      \           /                     |  
    %    |  |                       \_ _ _ _ _/                      |
    %    Ly |                        -L_Gapx-> ^                     |
    %    |  |                                  | L_Gapy              |
    %    |  |                <-Lix-> _ _ _ _ _ |                     | 
    %    |  |                ^    ^ /         \                      |     
    %    |  |                |   / /           \                     |
    %    |  |            Liy | Li /             \                    |
    %    |  |                | / /               \                   |
    %    |  |________________|/_/                 \__________________|
    %       <---------Lw-------><------ Lb -------><-----Lw--------->
    %     
    %                             Figure 2
    
    
    %%
    %{
    ***********************************************************************
        Parameters for 1st & 2nd Point Coordinates
    ***********************************************************************
    
    We calculate the parameters needed to find the coordinates of the first
    and second points of the elements on the channel.
    
    Polygon refers to the channel.
    %}
        
    %{
    Lenght of the first(also, the last) portion of the horizontal 
    bottom wall

    w: Wall
    %}
    mesh.L_Wall = (mesh.Lx - mesh.L_Base) / 2;

    %{
    Contraction Height

    i: Inclination
    %}
    mesh.L_Incx = (mesh.L_Base - mesh.L_Gapx) / 2;

    % Contraction Height
    mesh.L_Incy = (mesh.Ly - mesh.L_Gapy) / 2;

    % Coordinates of the vertices of the polygon
    mesh.X_Vc = [0                               0                    ;...
         mesh.L_Wall                             0                    ;...
         mesh.L_Wall + mesh.L_Incx               mesh.L_Incy          ;...
         mesh.L_Wall + mesh.L_Base - mesh.L_Incx mesh.L_Incy          ;...
         mesh.L_Wall + mesh.L_Base               0                    ;...
         2 * mesh.L_Wall  + mesh.L_Base          0                    ;...
         2 * mesh.L_Wall  + mesh.L_Base          mesh.Ly              ;...
         mesh.L_Wall + mesh.L_Base               mesh.Ly              ;...
         mesh.L_Wall + mesh.L_Base - mesh.L_Incx mesh.Ly - mesh.L_Gapy;...
         mesh.L_Wall + mesh.L_Incx               mesh.Ly - mesh.L_Gapy;...
         mesh.L_Wall                             mesh.Ly              ;...
         0                                       mesh.Ly              ];
    
    % Number of sides (line segments) of the polygon
    mesh.N_S_Poly = size(mesh.X_Vc, 1);
    
    %{
    Length of each side of polygon in a counter-clockwise direction
    
    S : Segment/Side
    
    L : Length
    %}
    mesh.S_L_Poly =...
        hypot(mesh.X_Vc([2 : end 1], 1) - mesh.X_Vc(:, 1),...
              mesh.X_Vc([2 : end 1], 2) - mesh.X_Vc(:, 2));
    
    %{
    Vector containing the number of elements at each line segment of
    the polygon
              
    Rounding up since we want at least one element on a side
    %}
    mesh.N_E_Poly = ceil(mesh.S_L_Poly / mesh.E_Size);
    
    %{
    Increasing the number of elements on the sensitive sides of the polygon
    %}
    mesh.N_E_Poly([2 3 9]) = 4 * mesh.N_E_Poly([2 3 9]);
    
    % Total number of elements on the channel wall
    mesh.N_E_0 = sum(mesh.N_E_Poly);
    
    %{
    Matrix containing the indices of the last elements of each line segment
    %}
    mesh.lastE_0 = cumsum(mesh.N_E_Poly);

    % Matrix containing the indices of the first elements of each circle
    mesh.firstE_0 = mesh.lastE_0 - mesh.N_E_Poly + 1;
    
    
    %%
    %{
    ***********************************************************************
        1st & 2nd Point Coordinates
    ***********************************************************************
    
    We calculate the coordinates of the first and second points of the 
    elements on the channel.
    %}
    [mesh.X1_0, mesh.X2_0] = BEM_StokesPartMeshPoly...
    (mesh.N_S_Poly, mesh.N_E_Poly, mesh.N_E_0, mesh.firstE_0,...
    mesh.lastE_0, mesh.X_Vc);
    
    
    %%
    %{
    ***********************************************************************
        Assembly Parameters
    ***********************************************************************
    
    These are the parameters needed for the assembly of the system 
    matrices.
    %}
    [mesh.Tg_0, mesh.E_L_0, mesh.E_N_0, mesh.XM_0]...
        = BEM_StokesPartAssemPars(mesh.X1_0, mesh.X2_0);
    
    
    %%
    %{
    ***********************************************************************
        Type and Value of the Boundary Condition (BC)
    ***********************************************************************
    
    Type Definition:
    
    1: Neumann (traction)
    
    0: Dirichlet (velocity)
    
    
    4 entries are reserved for each element in the matrix BC which is
    (2 * N) x 2.

    ie. 
        BC(1, 1) = 0  =>  x component of velocity will be specified at
                         node 1
        BC(1, 1) = 1  =>  x component of traction will be specified at
                         node 1
        BC(1, 2) = 0  =>  x component (of velocity or traction as 
                         specified in BC(1, 1)) is 0 at node 1
        BC(1, 2) = 1  =>  x component (of velocity or traction as 
                         specified in BC(1, 1)) is 1 at node 1
        BC(2, 1) = 0  =>  y component of velocity will be specified at
                         node 1
        BC(2, 1) = 1  =>  y component of traction will be specified at
                         node 1
        BC(2, 2) = 0  =>  y component (of velocity or traction 
                         depending on BC(2, 1)) is 0 at node 1
        BC(2, 2) = 1  =>  y component (of velocity or traction 
                         depending on BC(2, 1)) is 1 at node 1
    
    The above examples all indicate the boundary conditions for the first
    node. Starting from the BC(3, 3) entry, the other node specifications
    follow.
	%}
    
    %{
        Initializing the boundary condition matrix
    
    Since we fill up the matrix with zeros, this means all the elements
    have no slip and no penetration conditions on them.
    
    We then can change the inlet and outlet conditions only and the rest
    will stay the same
    %}
    mesh.BC_0 = zeros(2 * mesh.N_E_0, 2);
    
    %{
    Indicating that the x component of the velocity will be specified
    %}
	mesh.BC_0(2 * mesh.firstE_0(end) - 1 : 2 : 2 * mesh.lastE_0(end), 1)...
        = 0;
    
    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is vIn
    %}
    mesh.BC_0(2 * mesh.firstE_0(end) - 1 : 2 : 2 * mesh.lastE_0(end), 2)...
        = mesh.vIn;
    
    %{
    Indicating that the x component of the traction will be specified
    %}
	mesh.BC_0(2 * mesh.firstE_0(6) - 1 : 2 : 2 * mesh.lastE_0(6), 1) = 1;
    
    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is tOut
    %}
    mesh.BC_0(2 * mesh.firstE_0(6) - 1 : 2 : 2 * mesh.lastE_0(6), 2)...
        = mesh.tOut;
end