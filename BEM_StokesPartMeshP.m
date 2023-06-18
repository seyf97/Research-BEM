function mesh = BEM_StokesPartMeshP(mesh,values)
    % Particle Boundary Mesh Generation
    % Discretization of the circular boundaries
    % Forming the parameters related to the particle that are necessary for
    % the assembly of the system matrices


    % In case the particles are molecules
    if mesh.isMolUni
        mesh = BEM_StokesPartMeshMol(mesh);
    elseif mesh.isMolRand
        mesh = BEM_StokesPartMeshMolRand(mesh,values);
        
    % In case particles are rigid and circular
    else
        %%
        %{
        *******************************************************************
            Coordinates of the first points of the elements

        P : abbreviation for Particle

        X1_P : N_Ep x 2 matrix
        *******************************************************************
        %}
        
        %{
        Initializing the matrix storing the first points of the elements of
        the circles
        %}
% % % % % % % % % % % % % % % % % % % % %         mesh.X1_P = zeros(mesh.N_E_P_O, 2);
% % % % % % % % % % % % % % % % % % % % %         
% % % % % % % % % % % % % % % % % % % % %         %{
% % % % % % % % % % % % % % % % % % % % %         Initializing the matrix for the 2nd points the same way as the 1st 
% % % % % % % % % % % % % % % % % % % % %         one
% % % % % % % % % % % % % % % % % % % % %         %}
% % % % % % % % % % % % % % % % % % % % %         mesh.X2_P = zeros(mesh.N_E_P_O, 2);

        %
            %
            [tempX1, tempX2] = BEM_StokesPartMeshCircle...
             (mesh.Rad, mesh.X_CG, mesh.rotAng, mesh);

            %
            mesh.X1_P(mesh.firstE_P : mesh.lastE_P, :) = tempX1;

            %
            mesh.X2_P(mesh.firstE_P : mesh.lastE_P, :) = tempX2;
    end


    %%
    %{
    *******************************************************************
        Particle Boundary BC
    *******************************************************************
    %}

    % Initializing the circular boundary condition matrix
    mesh.BC_P = zeros(2 * sum(mesh.N_E_P), 2);

    %{
    Indicating that the x component of the velocity will be specified
    %}
    mesh.BC_P(1 : 2 : end, 1) = 0;

    %{
    Indicating that the x component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    mesh.BC_P(1 : 2 : end, 2) = 0;

    %{
    Indicating that the y component of the velocity will be specified
    %}
    mesh.BC_P(2 : 2 : end, 1) = 0;

    %{
    Indicating that the y component (of velocity or traction 
    depending on the boundary condition type specified above) is 0
    %}
    mesh.BC_P(2 : 2 : end, 2) = 0;


    %%
    %{
    *******************************************************************
        Assembly Parameters
    *******************************************************************

    These are the parameters needed for the assembly of the system 
    matrices.
    %}
    [mesh.Tg_P, mesh.E_L_P, mesh.E_N_P, mesh.XM_P]...
        = BEM_StokesPartAssemPars(mesh.X1_P, mesh.X2_P);


    %%
    %{
    *******************************************************************
        Relative Position Matrix of Element Nodes with respect to CG
    *******************************************************************

    Finding the differences between the coordinates of the midpoints 
    of the elements on the particle and the coordinates of the center
    of gravity of the paticle

    S : abbreviation for Segment

    Tg : abbreviation for Tangent

    STg : N_Ep x 2 matrix
    %}

    % Initializing the relative position matrix
    mesh.S_Tg = zeros(sum(mesh.N_E_P), 2);

    %{
    Using the same index matrix used for filling up the coordinates 
    matrices
    %}
    for i = 1 : mesh.N_P

        mesh.S_Tg(mesh.lastE_P(i) - mesh.N_E_P(i) + 1 :...
            mesh.lastE_P(i), :) =...
        mesh.XM_P(mesh.lastE_P(i) - mesh.N_E_P(i) + 1 :...
            mesh.lastE_P(i), :) - mesh.X_CG(i, :);
    end


    %%
    %{
    *******************************************************************
        Line Segment Length
    *******************************************************************

    Distance (length of the line segment) from the computational nodes
    of the elements on the particle surface to the CG of the particle

    We find it by using the Pythagorean Theorem.

    SL : abbreviation for Line Segment from the CG to the arbitrary 
         node

    S : N_Ep x 1
    *******************************************************************
    %}

    mesh.S_L = hypot(mesh.S_Tg(:, 1), mesh.S_Tg(:, 2));


    %%
    %{
    *******************************************************************
        Unit Vector Normal to Line Segment
    *******************************************************************

    Finding the components of the unit vector that is normal to the 
    line segment drawn from the CG of the particle to the node (the 
    direction of the normal is selected using right hand rule, 
    counter-clockwise).

    N : abbreviation for the aforementioned normal unit vector
    %}

    % x coordinates
    mesh.S_N(:, 1) = - mesh.S_Tg(:, 2) ./ mesh.S_L;

    % y coordinates
    mesh.S_N(:, 2) = mesh.S_Tg(:, 1) ./ mesh.S_L;
    
    
    %%
    %{
    ***********************************************************************
    Addition of Particle Parameters to Channel Parameters
    ***********************************************************************
    %}

    %
    mesh.X1 = [mesh.X1_0; mesh.X1_P];

    %
    mesh.X2 = [mesh.X2_0; mesh.X2_P];

    %
    mesh.XM = [mesh.XM_0; mesh.XM_P];

    %
    mesh.E_N = [mesh.E_N_0; mesh.E_N_P];

    %
    mesh.E_L = [mesh.E_L_0; mesh.E_L_P];

    %
    mesh.Tg = [mesh.Tg_0; mesh.Tg_P];
    
    %
    mesh.BC = [mesh.BC_0; mesh.BC_P];


end