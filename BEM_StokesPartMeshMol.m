function mesh = BEM_StokesPartMeshMol(mesh) 
    % Boundary mesh for attached particles
    %
    
    %CIRCULAR ELEMENTS
    
    %{
    Angle formed by the sector that spans the main particle in between
    the two points at which an arbitrary small particle is attached to it
    
    Note that this angle will be the same for all the attached particles.
    %}
    mesh.angSmall = 2 * acos((mesh.Rad ^ 2 + (mesh.Rad + mesh.RadSmall -...
        mesh.penAm) ^ 2 - mesh.RadSmall ^ 2) / (2 * mesh.Rad * ...
        (mesh.Rad + mesh.RadSmall - mesh.penAm)));
    
    %Lambda=Angle btw the centers of two small circles
    angCG_Small=2*pi/mesh.N_P_Small;

    %Creating the vector to store the centers of small circles
    mesh.X_CG_Small=zeros(mesh.N_P_Small,2);

    %Half the distance between the two coincident points of the Small cirlce
    %(S) and the Large circle (L)
    b=mesh.Rad*real(sin(mesh.angSmall));

    %Distance between the center of a small circle and Origin
    L=mesh.Rad*real(cos(mesh.angSmall)) + sqrt(mesh.RadSmall.^2 - b.^2);
    
    %Dist btw outer points of the L and S circle
    SOuterPt=L-mesh.RadSmall;
    LOuterPt=mesh.Rad;
    mesh.Dist=LOuterPt-SOuterPt;


    %Creating the X coordinates of the CentofSmall circles. Afterwards,
    %deleting the last element since it overlaps the first
    dummy = L*cos(mesh.rotAng+0:-2*pi/mesh.N_P_Small:mesh.rotAng-2*pi);
    dummy(end)=[];
    
    %Putting the X coordinates of the Center of Small circles into a matrix
    mesh.X_CG_Small(:,1) = dummy;
    
    %Creating the Y coordinates of the CentofSmall circles. Again the last
    %element will be deleted
    dummy=L*sin(mesh.rotAng+0:-2*pi/mesh.N_P_Small:mesh.rotAng-2*pi);
    dummy(end)=[];
    
    %Putting the Y coordinates into the matrix
    mesh.X_CG_Small(:,2)=dummy;
    
    %Beta, the angle of the starting point of meshing of the small circle
    Beta=atan(sqrt(mesh.RadSmall.^2 - b.^2)/b);
    
    %Initializing the Matrix that will keep the small Circle coordinates
    %(INCLUDING THE TWO INTERCEPTION POINTS)
    XSmallCircle = zeros((mesh.N_E_P_Small+1)*mesh.N_P_Small,2);
    
    %Number of elements on a single gap of the Large Circle between two small circles
    GapNumEl=mesh.N_E_P_Large/mesh.N_P_Small;

    %Initializing the coordinates of the Large Circle. 
    %(INCLUDING THE TWO INTERCEPTION POINTS)
    XLargeCircle=zeros((GapNumEl+1)*mesh.N_P_Small,2);

    %The angle of a single gap which will be meshed
    MeshAngle=angCG_Small-2*mesh.angSmall;
    %Angle of increment for the next gap
    IncrAngle=MeshAngle+2*mesh.angSmall;

    %Initialization of all Coordinate Points, advancing in clockwise
    AllPoints=zeros(mesh.N_P_Small*(GapNumEl+mesh.N_E_P_Small),2);



    for i=1:mesh.N_P_Small
            %            --------------- 1) SMALL CIRCLE -----------------------
        %First we create the small circle coordinates
        %Adding Lambda to the angle (clockwise) in order to change the circle
        %in every iteration

        x=mesh.RadSmall*cos(linspace(mesh.rotAng+(pi/2)+Beta-angCG_Small*(i-1),mesh.rotAng+(-pi/2)-Beta-angCG_Small*(i-1),mesh.N_E_P_Small+1));
        y=mesh.RadSmall*sin(linspace(mesh.rotAng+(pi/2)+Beta-angCG_Small*(i-1),mesh.rotAng+(-pi/2)-Beta-angCG_Small*(i-1),mesh.N_E_P_Small+1));

        %Moving the small circles to their center of gravities
        x=x+mesh.X_CG_Small(i,1);
        y=y+mesh.X_CG_Small(i,2);

        %Making sure that all coordinates are consecutive
        XSmallCircle(((i-1)*(mesh.N_E_P_Small+1)+1):(mesh.N_E_P_Small+1)*i,1)=x;
        XSmallCircle(((i-1)*(mesh.N_E_P_Small+1)+1):(mesh.N_E_P_Small+1)*i,2)=y;

        %Taking out the interception points with the L circle
        x(1)=[];
        x(end)=[];
        y(1)=[];
        y(end)=[];

        %Adding all the points to the main matrix
        AllPoints(((1+(GapNumEl+mesh.N_E_P_Small)*(i-1)):(GapNumEl+mesh.N_E_P_Small)*(i-1)+(mesh.N_E_P_Small-1)),1) = x;
        AllPoints(((1+(GapNumEl+mesh.N_E_P_Small)*(i-1)):(GapNumEl+mesh.N_E_P_Small)*(i-1)+(mesh.N_E_P_Small-1)),2) = y;





        %            --------------- 2) LARGE CIRCLE -----------------------

        %Now we create the Large Circle coordinates
        %Adding Incr angle for every consecutive gaps (clockwise)

        x=mesh.Rad*cos(linspace(mesh.rotAng-mesh.angSmall-IncrAngle*(i-1),mesh.rotAng-mesh.angSmall-MeshAngle-IncrAngle*(i-1),GapNumEl+1));
        y=mesh.Rad*sin(linspace(mesh.rotAng-mesh.angSmall-IncrAngle*(i-1),mesh.rotAng-mesh.angSmall-MeshAngle-IncrAngle*(i-1),GapNumEl+1));

        %Making sure that all coordinates are consecutive
        %(INCLUDING THE TWO INTERCEPTION POINTS)
        %
        XLargeCircle(((i-1)*(GapNumEl+1)+1):(GapNumEl+1)*i,1)=x;
        XLargeCircle(((i-1)*(GapNumEl+1)+1):(GapNumEl+1)*i,2)=y;

        %Adding the Gap elements in the all points matrix
        AllPoints((mesh.N_E_P_Small)+(i-1)*(mesh.N_E_P_Small+GapNumEl):(mesh.N_E_P_Small+GapNumEl)+(i-1)*(mesh.N_E_P_Small+GapNumEl),1)=x;
        AllPoints((mesh.N_E_P_Small)+(i-1)*(mesh.N_E_P_Small+GapNumEl):(mesh.N_E_P_Small+GapNumEl)+(i-1)*(mesh.N_E_P_Small+GapNumEl),2)=y;

    end

    %Creating the X1 Points of the Particle Elements
    X1=AllPoints;
    X1=X1+repmat(mesh.X_CG,size(X1,1),1);

    %Creating the X2 Points of the Particle Elements
    X2=zeros(size(AllPoints,1),2);
    X2(2:end,1)=X1(1:(end-1),1);
    X2(1,1)=X1(end,1);
    X2(2:end,2)=X1(1:(end-1),2);
    X2(1,2)=X1(end,2);

    %DUZELTME
    dummy=X1;
    X1=X2;
    X2=dummy;

    mesh.X1_P = X1;

    mesh.X2_P = X2;
    
    mesh.X_CG_Small = mesh.X_CG_Small + mesh.X_CG;
    
    
    %%
    %{
    ***********************************************************************
        Assembly Parameters
    ***********************************************************************
 
    These are the parameters needed for the assembly of the system 
    matrices.
    %}
    [mesh.Tg_P, mesh.E_L_P, mesh.E_N_P, mesh.XM_P]...
        = BEM_StokesPartAssemPars(mesh.X1_P, mesh.X2_P);
    
    
    %%
    %{
    ***********************************************************************
        Relative Position Matrix of Element Nodes with respect to CG
    ***********************************************************************
    
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
        
        mesh.S_Tg(mesh.lastE_P(i) - mesh.N_E_P(i) + 1 : mesh.lastE_P(i), :) =...
            mesh.XM_P(mesh.lastE_P(i) - mesh.N_E_P(i) + 1 : mesh.lastE_P(i), :) - mesh.X_CG(i, :);
    end
    
    
    %%
    %{
    ***********************************************************************
        Line Segment Length
    ***********************************************************************
    
    Distance (length of the line segment) from the computational nodes
    of the elements on the particle surface to the CG of the particle
    
    We find it by using the Pythagorean Theorem.
    
    SL : abbreviation for Line Segment from the CG to the arbitrary node
    
    S : N_Ep x 1
    ***********************************************************************
    %}
    
    mesh.S_L = hypot(mesh.S_Tg(:, 1), mesh.S_Tg(:, 2));
    
    
    %%
    %{
    ***********************************************************************
        Unit Vector Normal to Line Segment
    ***********************************************************************
    
    Finding the components of the unit vector that is normal to the line 
    segment drawn from the CG of the particle to the node (the direction of
    the normal is selected using right hand rule, counter-clockwise).
    
    N : abbreviation for the aforementioned normal unit vector
    %}
    
    % x coordinates
    mesh.S_N(:, 1) = - mesh.S_Tg(:, 2) ./ mesh.S_L;
    
    % y coordinates
    mesh.S_N(:, 2) = mesh.S_Tg(:, 1) ./ mesh.S_L;
    
end