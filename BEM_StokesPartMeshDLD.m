function mesh = BEM_StokesPartMeshDLD(mesh)
    %
    
    % Number of elements on a single post
    mesh.N_E_Post = 3 * ceil(2 * pi * mesh.postRad / mesh.E_Size);

    %
    [mesh.X1_Post, mesh.X2_Post, mesh.X_CG_Post,...
        mesh.N_Post, mesh.N_E_Post_O] =...
        BEM_StokesPartMeshPostArr...
        (mesh.offsetx, mesh.offsety, mesh.postRad, mesh.shifty,...
        mesh.postArrLx, mesh.postArrLy, mesh.X_CG_1, mesh.N_E_Post);
    
    
    %%
    %{
    ***********************************************************************
        Assembly Parameters
    ***********************************************************************
    
    These are the parameters needed for the assembly of the system 
    matrices.
    %}
    [mesh.Tg_Post, mesh.E_L_Post, mesh.E_N_Post, mesh.XM_Post]...
        = BEM_StokesPartAssemPars(mesh.X1_Post, mesh.X2_Post);
    
    
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
    %}
    mesh.BC_Post = zeros(2 * mesh.N_E_Post_O, 2);
    
    
    
    %%
    %{
    ***********************************************************************
        Addition of Posts
    ***********************************************************************
    
    Adding the post parameters needed for the assembly of the system 
    matrices
    %}
    
    % 1st Points Matrix
    mesh.X1_0 = [mesh.X1_0; mesh.X1_Post];
    
    % 2nd Points Matrix
    mesh.X2_0 = [mesh.X2_0; mesh.X2_Post];
    
    % Middle Points (Computational Nodes) Matrix
    mesh.XM_0 = [mesh.XM_0; mesh.XM_Post];
    
    % Element Lengths Vector
    mesh.E_L_0 = [mesh.E_L_0; mesh.E_L_Post];
    
    % Normal Vectors Matrix
    mesh.E_N_0 = [mesh.E_N_0; mesh.E_N_Post];
    
    % Relative Position Matrix
    mesh.Tg_0 = [mesh.Tg_0; mesh.Tg_Post];
    
    % Boundary Conditions Matrix
    mesh.BC_0 = [mesh.BC_0; mesh.BC_Post];
    
    % Number of Elements on the channel including the posts
    mesh.N_E_0 = mesh.N_E_0 + mesh.N_E_Post_O;
end