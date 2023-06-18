function mesh = BEM_StokesPartMeshPFF(mesh)
    % Discretization of the channel wall
    %       
    % 
    %                                    _________________  
    %              |\                   /                 |  ^
    %              | \                 /                  |  |
    %              \  \               /                   |  |
    %               \  \             /                    |  |
    %                \  \           /)                    |  |
    %                 \  \         /   )                  |  |
    %                  \  \_______/      )                |  |
    %                 ( \  -L_Pinx-> ^     )              |  |
    %          InAng (   \           |      ) bounAng     |  L_Broy
    %                (   /    L_Piny |      )             |  |
    %                 ( /  _______   |     )              |  |
    %                  /  /       \      )                |  |  
    %                 /  /         \   )                  |  |
    %                /  /           \)                    |  |
    %               /  /             \                    |  |
    %              /  /               \                   |  |
    %            ^ | /                 \                  |  |
    %      L_Iny | |/                   \_________________|  |
    %                                    <---- L_Brox ---->
    %     
    %                             Figure 3
    
    
    %%
    %{
    ***********************************************************************
        Parameters for 1st & 2nd Point Coordinates
    ***********************************************************************
    
    We calculate the parameters needed to find the coordinates of the first
    and second points of the elements on the channel.
    
    Polygon refers to the channel.
    %}
    
    %
    mesh.L_PFF_1x = (mesh.L_Broy - mesh.L_Piny) /...
        (2 * tan(mesh.nozzleAng / 2));
    
    %
    mesh.L_PFF_1y = (mesh.L_Broy - mesh.L_Piny) / 2;
    
    %
    mesh.L_PFF_2x = mesh.L_PFF_1x + mesh.L_Pinx;

    %
    mesh.L_PFF_3x = mesh.L_PFF_2x + mesh.L_Broy /...
        (2 * tan(mesh.diffuserAng / 2));

    %
    mesh.L_PFF_4x = mesh.L_PFF_3x + mesh.L_Brox;

    %
    mesh.L_PFF_5y = (mesh.L_Broy + mesh.L_Piny) / 2;

    %
    mesh.L_PFF_6y = mesh.L_Broy - mesh.L_Iny;

    %
    mesh.L_PFF_5x = (mesh.L_Broy / 2 - mesh.L_Iny) /...
        tan(mesh.nozzleAng / 2);

    % Coordinates of the vertices of the polygon
    mesh.X_Vc = [0                   0              ;...
                 mesh.L_PFF_1x       mesh.L_PFF_1y  ;...
                 mesh.L_PFF_2x       mesh.L_PFF_1y  ;...
                 mesh.L_PFF_3x       0              ;...
                 mesh.L_PFF_4x       0              ;...
                 mesh.L_PFF_4x       mesh.L_Broy    ;...
                 mesh.L_PFF_3x       mesh.L_Broy    ;...
                 mesh.L_PFF_2x       mesh.L_PFF_5y  ;...
                 mesh.L_PFF_1x       mesh.L_PFF_5y  ;...
                 0                   mesh.L_Broy    ;...
                 0                   mesh.L_PFF_6y  ;...
                 mesh.L_PFF_5x       mesh.L_Broy / 2;...
                 0                   mesh.L_Iny     ];
    
    % Number of sides (line segments) the polygon
    mesh.N_S_Poly = size(mesh.X_Vc, 1);
    
    %{
    Length of each side of polygon except for the last side in a
    counter-clockwise direction
    
    S : Segment/Side
    
    L : Length
    %}
    mesh.S_L_Poly =...
        hypot(mesh.X_Vc([2 : end 1], 1) - mesh.X_Vc(:, 1),...
              mesh.X_Vc([2 : end 1], 2) - mesh.X_Vc(:, 2));
    
    %{
    Vector containing the number of elements of each line segment of
    the polygon
    %}
    mesh.N_E_Poly = ceil(mesh.S_L_Poly / mesh.E_Size);
    
    % Increasing the number of elements on the sensitive sides
    mesh.N_E_Poly([7 8]) = 1 * [mesh.N_E_Poly(7) mesh.N_E_Poly(8)];
    
    % Total number of elements on the channel
    mesh.N_E_0 = sum(mesh.N_E_Poly);
    
    %{
    Matrix containing the indices of the last elements of each line segment
    %}
    mesh.lastE_0 = cumsum(mesh.N_E_Poly);
    
    %{
    Matrix containing the indices of the first elements of each side of
    the polygon
    %}
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
    Indicating the x component (of velocity or traction 
    depending on the boundary condition type specified above)
    %}
    mesh.BC_0(2 * mesh.firstE_0(end) - 1 : 2 : 2 * mesh.lastE_0(end), 2)...
        = mesh.vIn0 * cos(mesh.nozzleAng / 2);
    
    %{
    Indicating that the y component of the velocity will be specified
    %}
	mesh.BC_0(2 * mesh.firstE_0(end) : 2 : 2 * mesh.lastE_0(end), 1) = 0;
    
    %{
    Indicating the y component (of velocity or traction 
    depending on the boundary condition type specified above)
    %}
    mesh.BC_0(2 * mesh.firstE_0(end) : 2 : 2 * mesh.lastE_0(end), 2) =...
        mesh.vIn0 * sin(mesh.nozzleAng / 2);
    
    %{
    Indicating that the x component of the velocity will be specified
    %}
	mesh.BC_0(2 * mesh.firstE_0(end - 3) - 1 : 2 :...
              2 * mesh.lastE_0(end - 3), 1) = 0;
    
    %{
    Indicating the x component (of velocity or traction 
    depending on the boundary condition type specified above)
    %}
    mesh.BC_0(2 * mesh.firstE_0(end - 3) - 1 : 2 :...
       2 * mesh.lastE_0(end - 3), 2) = mesh.vInP * cos(mesh.nozzleAng / 2);
    
    %{
    Indicating that the y component of the velocity will be specified
    %}
	mesh.BC_0...
     (2 * mesh.firstE_0(end - 3) : 2 : 2 * mesh.lastE_0(end - 3), 1) = 0;
    
    %{
    Indicating that the y component (of velocity or traction 
    depending on the boundary condition type specified above) is vIn
    %}
    mesh.BC_0...
     (2 * mesh.firstE_0(end - 3) : 2 : 2 * mesh.lastE_0(end - 3), 2) =...
     - mesh.vInP * sin(mesh.nozzleAng / 2);
    
    %{
    Indicating that the x component of the traction will be specified
    %}
	mesh.BC_0(2 * mesh.firstE_0(5) - 1 : 2 : 2 * mesh.lastE_0(5), 1) = 1;
    
    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is tOut
    %}
    mesh.BC_0(2 * mesh.firstE_0(5) - 1 : 2 : 2 * mesh.lastE_0(5), 2) =...
     mesh.tOut;
end