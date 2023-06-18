function mesh = BEM_StokesPartMeshSRC(mesh)
    % Discretization of the channel wall
    % 
    % 
    % y axis ^
    %        ^
    %        |  <------------------------   Lx  ---------------------->  
    %        |   ______________________________________________________
    %        ^  |             ^                                        | 
    %        |  |   <   >     |                                        | 
    %        |  |  <  .  >    | b                                      | 
    %        |  |   <   >     |                                        | 
    %        |  |             |                                        |
    %        |  |    <   >    |                                        | 
    %        |  |  <       >  |                                        | 
    %        |  | <   .- a -> |                                        | 
    %        |  |  <       >  |                                        | 
    %        |  |    <   >    |                                        | 
    %        |  |             |                                        |
    %        |  |_._._._._._._|._._._._._._._._._._._._._._._._._._._._|
    %        |  |                                    mid line          |
    %        |  |    <   >                                             | 
    %        |  |  <       >                                           | 
    %        |  | <   .- a ->                                          | 
    %        |  |  <       >                                           | 
    %        |  |    <   >                                             |
    %        |  |                                                      |
    %        |  |    <   >                                             | 
    %        |  |   <  .  >                                            | 
    %        |  |    <   >                                             | 
    % Origin ._ |______________________________________________________|->> 
    %                                                                x axis
    %                         
    %                             Figure 1

    
    
    %%
    %{
    ***********************************************************************
        1st & 2nd Point Coordinates
    ***********************************************************************
    
    We calculate the coordinates of the first and second points of the 
    elements on the channel.
    
    Polygon refers to the channel in this function.
    %}
    
    %{
        Number of elements on the bottom or top sides of the channel wall

    If there is a contraction (for performing hydrodynamic separation)in 
    the channel, the number of elements here should be at least 5. 
    Otherwise, it should be at least 1.

    N : Number

    E : Element

    0 : Channel Wall

    x : horizontal dimension (as indicated in the Figure 1 inside the Mesh0
        function)

    N_E0x : 1 x 1
    %}
    mesh.N_E_0x = ceil(mesh.Lx / mesh.E_Size);
    
    %{
        Number of elements on the left or right sides of the channel wall

    If there is a contraction (for performing hydrodynamic separation)in 
    the channel, the number of elements here should be at least 5. 
    Otherwise, it should be at least 1.

    y : vertical dimension (as indicated in the Figure 1 under the section 
        Geometry Specifications)

    N_E0y : 1 x 1
    %}
    mesh.N_E_0y = ceil(mesh.Ly / mesh.E_Size);
    
    %{
    Vector containing the number of elements of each line segment of
    the polygon
    %}
    mesh.N_E_Poly = [mesh.N_E_0x;...
                     mesh.N_E_0y;...
                     mesh.N_E_0x;...
                     mesh.N_E_0y];
    
    % Total number of elements on the channel
    mesh.N_E_0 = sum(mesh.N_E_Poly);
    
    % Coordinates of the vertices of the polygon
    mesh.X_Vc =      [0 0      ;...
                mesh.Lx 0      ;...
                mesh.Lx mesh.Ly;...
                      0 mesh.Ly];
    
    % Number of sides (line segments) the polygon
    mesh.N_S_Poly = size(mesh.X_Vc, 1);
    
    %{
        Matrix containing the indices of the last elements of each line
        segment
    %}
    mesh.lastE_0 = cumsum(mesh.N_E_Poly);
    
    % Matrix containing the indices of the first elements of each circle
    mesh.firstE_0 = mesh.lastE_0 - mesh.N_E_Poly + 1;
    
    % Forming the first and second points of the elements of the boundary
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
    -----------------------------------------------------------------------
        Bottom boundary
    %}
    
    % Initializing the bottom boundary condition matrix
    BC_B = zeros(2 * mesh.N_E_0x, 2);
    
    %{
    Indicating that the x component of the velocity will be specified
    %}
    BC_B(1 : 2 : end, 1) = 0;
    
    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    BC_B(1 : 2 : end, 2) = 0;
    
    %{
    Indicating that the y component of the velocity will be specified
    %}
    BC_B(2 : 2 : end, 1) = 0;
    
    %{
    Indicating that the y component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    BC_B(2 : 2 : end, 2) = 0;
    %----------------------------------------------------------------------
    
    %{
    -----------------------------------------------------------------------
        Right boundary
    %}
    
    % Initializing the right boundary condition matrix
    BC_R = zeros(2 * mesh.N_E_0y, 2);
    
    %{
    Indicating that the x component of the traction will be specified
    %}
	BC_R(1 : 2 : end, 1) = 1;
    
    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    BC_R(1 : 2 : end, 2) = mesh.tOut;
    
    %{
    Indicating that the y component of the velocity will be specified
    %}
    BC_R(2 : 2 : end, 1) = 0;
    
    %{
    Indicating that the y component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    BC_R(2 : 2 : end, 2) = 0;
    %----------------------------------------------------------------------
    
    %{
    -----------------------------------------------------------------------
        Top boundary
    
    For illustration purposes we impose the top boundary condition type and
    value as well though we could just as easily made it equal to the
    specifications about the bottom boundary in a symmetric configuration
    to save from the computational cost
    %}
    
    % Initializing the top boundary condition matrix
    BC_T = zeros(2 * mesh.N_E_0x, 2);
    
    %{
    Indicating that the x component of the velocity will be specified
    %}
    BC_T(1 : 2 : end, 1) = 0;
    
    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    BC_T(1 : 2 : end, 2) = 0;
    
    %{
    Indicating that the y component of the traction will be specified
    %}
    BC_T(2 : 2 : end, 1) = 0;
    
    %{
    Indicating that the y component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    BC_T(2 : 2 : end, 2) = 0;
    %----------------------------------------------------------------------
    
    %{
    -----------------------------------------------------------------------
        Left boundary
    %}
    
    % Initializing the left boundary condition matrix
    BC_L = zeros(2 * mesh.N_E_0y, 2);
    
    %{
    Indicating that the x component of the velocity will be specified
    %}
	BC_L(1 : 2 : end, 1) = 0;
    
    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is vIn
    %}
    BC_L(1 : 2 : end, 2) = mesh.vIn;
    
    %{
    Indicating that the y component of the velocity will be specified
    %}
    BC_L(2 : 2 : end, 1) = 0;
    
    %{
    Indicating that the y component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    BC_L(2 : 2 : end, 2) = 0;
    %----------------------------------------------------------------------
    
    %{
    Creating the matrix which specifies the type of the boundary condition
    at each element in the counter-clockwise direction along the sides of
    the rectangle starting from the origin and going in the clockwise 
    direction along the circle inside
    %}
	mesh.BC_0 = [BC_B; BC_R; BC_T; BC_L];
end