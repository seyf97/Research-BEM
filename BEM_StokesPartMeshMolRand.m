function mesh = BEM_StokesPartMeshMolRand(mesh,values) 

%SEYFO 20/10/18
%01:58


theta = mesh.rotAng;
elSize = mesh.E_Size/mesh.E_C_Multiplier;

    %{
    Angle formed by the sector that spans the main particle in between
    the two points at which an arbitrary small particle is attached to it
    
    Note that this angle will be the same for all the attached particles.
    %}

mesh.angSmall = 2 * acos((mesh.Rad ^ 2 + (mesh.Rad + mesh.RadSmall -...
        mesh.penAm) ^ 2 - mesh.RadSmall ^ 2) / (2 * mesh.Rad * ...
        (mesh.Rad + mesh.RadSmall - mesh.penAm)));
    
    %Converting mesh.angSmall into degrees since i felt like using degrees
    %in the following parts
    alpha = mesh.angSmall*(180/pi);
    
numOfMolecules = mesh.N_P_Small;
%The minimum angle needed in between two circles in order to prevent
%interfering.
minAngle = asin( mesh.RadSmall / (mesh.Rad+mesh.RadSmall) );
minAngle = 2*ceil ( (180/pi) * minAngle );

dummyVec = zeros(1,numOfMolecules);

% VER: 1 as of 16.10 15:22
% VER: 2 as of 20.10 00:30


%Statement to make sure that the mesh is created only once because its
%random
%When the function is called next time, it will directly skip creating the
%points and will only translate and then rotate the points by taking the
%realCG as the center.

if ~mesh.iter
    
    %Counter to prevent the random loop to go forever
    cnt1 = 0;
    limitcnt1 = 5e6;
    
for i = 1:numOfMolecules
    
    dummyVec(i) = randi([0 359]);
    
    while (i ~= 1) && ( max(( abs(dummyVec(i) - dummyVec(i-1:-1:1)) <= minAngle )) || ...
            max(( abs(dummyVec(i) - dummyVec(i-1:-1:1)) >= (360 - minAngle) )) )
        
        dummyVec(i) = randi([0 359]);
        cnt1 = cnt1 + 1;
        %making sure that if the random gen can not fit the molecules, the
        %code gives an error message
        
        if cnt1 >  limitcnt1
            error('Max number of iterations (5e6) reached. Please restart the code.(Use less num of Molecules)');
            break
        end
        
            
    end

end


dummyVec = sort(dummyVec);
dummyVec = -dummyVec;
mesh.angCircles = dummyVec;










%Creating the vector to store the centers of small circles
CentSmallCircles=zeros(numOfMolecules,2);

%Half the distance between the two coincident points of the Small cirlce
%(S) and the Large circle (L)
b=mesh.Rad*sind(alpha);

%Distance between the center of a small circle and Origin
L=mesh.Rad*cosd(alpha) + sqrt(mesh.RadSmall.^2 - b.^2);

%Dist btw outer points of the L and S circle. Must be larger than 0.
SOuterPt=L-mesh.RadSmall;
LOuterPt=mesh.Rad;
mesh.dist=LOuterPt-SOuterPt;

    
        

for i = 1: numOfMolecules
    
    CentSmallCircles(i,1) = L*cosd(dummyVec(i));
    
    CentSmallCircles(i,2) = L*sind(dummyVec(i));
    
end

%Creating the matrix which will store the center of small circles
mesh.X_CG_Small = CentSmallCircles;

%Beta, the angle of the starting point of meshing of the small circle
Beta=atand(sqrt(mesh.RadSmall.^2 - b.^2)/b);
mesh.Beta = Beta;
% ( elSize/(2*mesh.RadSmall) )
% (2*asin( elSize/(2*mesh.RadSmall) ) )
% (pi + 2*Beta*(pi/180))

%Obtaining the number of elements on a small circle
numElS = ceil( (pi + 2*Beta*(pi/180)) / (2*asin( elSize/(2*mesh.RadSmall) ) ) );

% numElAllS = ceil( (2*pi) / (2*asin( elSize/(2*mesh.RadSmall) ) ) );

%Initializing the Matrix that will keep the small Circle coordinates
%(INCLUDING THE TWO INTERCEPTION POINTS)
XSmallCircle = zeros((numElS+1)*numOfMolecules,2);



%Finding the no. of elements on the gaps for every gap
GapEl = zeros(1,numOfMolecules);




%Initializing the vector which will keep the Lambda angles
lambdaVector = zeros(1,numOfMolecules);

for i = 1:numOfMolecules
    
    if i ~= numOfMolecules
        lambdaVector(i) = -dummyVec(i+1) + dummyVec(i);
        GapEl(i) = ceil( (lambdaVector(i)*(pi/180) - 2*alpha*(pi/180)) / (2*asin( elSize/(2*mesh.Rad) ) ) );
        
    else
       lambdaVector(i) = 360 + dummyVec(i) -dummyVec(1);
       GapEl(i) = ceil( (lambdaVector(i)*(pi/180) - 2*alpha*(pi/180)) / (2*asin( elSize/(2*mesh.Rad) ) ) );
          
    end

end

%NUMBER OF ELEMENTS
mesh.N_E_P_Large = sum(GapEl);
mesh.N_E_P_Small = numElS;
mesh.N_E_P_O = mesh.N_E_P_Large + mesh.N_E_P_Small*mesh.N_P_Small;
mesh.N_E_P = mesh.N_E_P_O;

%Matrix containing the indices of the last elements of each circle
mesh.lastE_P = mesh.N_E_P;

% Matrix containing the indices of the first elements of each circle
mesh.firstE_P = mesh.lastE_P - mesh.N_E_P + 1;


%Removing 2 points from gaps since the small circles include it
GapPts = GapEl - 2;
XGapPts = zeros(sum((GapPts)),2);

AllXPoints = [XSmallCircle ; XGapPts];
%Creating the coordinate points of the small circles and gaps

%The counter used when adding coordinates 
coordinateSum = 0;

for i = 1:numOfMolecules

    
   %Small Circles
   x=mesh.RadSmall*cosd(linspace((180/2)+Beta+dummyVec(i),(-180/2)-Beta+dummyVec(i),numElS+1));
   y=mesh.RadSmall*sind(linspace((180/2)+Beta+dummyVec(i),(-180/2)-Beta+dummyVec(i),numElS+1));
   
    x = x + CentSmallCircles(i,1);
    y = y + CentSmallCircles(i,2);
    
    x(end) = [];
    y(end) = [];
   
   
   %Large Circles
   
   X=mesh.Rad*cosd(linspace( dummyVec(i) - alpha , dummyVec(i) - lambdaVector(i) + alpha , GapEl(i) + 1));
   Y=mesh.Rad*sind(linspace( dummyVec(i) - alpha , dummyVec(i) - lambdaVector(i) + alpha , GapEl(i) + 1));
   
   X(end) = [];
   Y(end) = [];
   
   AllXPoints(coordinateSum + 1 : coordinateSum + (numElS + 1) + (GapEl(i) - 1),1) = [x X];
   AllXPoints(coordinateSum + 1 : coordinateSum + (numElS + 1) + (GapEl(i) - 1),2) = [y Y];    
   
   coordinateSum = coordinateSum + (numElS + 1) + (GapEl(i) - 1);
   
end




%FINDING THE FIRST REAL CENTER OF GRAVITY DUE TO THE RANDOM DIST OF
%MOLECULES

areaSingleCircle = pi*mesh.RadSmall.^2;
areaLargeCircle = pi*mesh.Rad.^2;

CentOfGravityX = ( areaSingleCircle*sum((CentSmallCircles(:,1))) )/( areaSingleCircle*numOfMolecules + areaLargeCircle);
CentOfGravityY = ( areaSingleCircle*sum((CentSmallCircles(:,2))) )/( areaSingleCircle*numOfMolecules + areaLargeCircle);

CentOfGrav = [CentOfGravityX CentOfGravityY];
x = CentOfGrav(1);
y = CentOfGrav(2);


%CREATING X1 X2 PTS AND CREATING THE REAL_CG. THEN TRANSLATING THEM TO
%THEIR INITIAL POSITION.

%Creating the X1 Points of the Particle Elements
X1 = AllXPoints;
%Translating the Particle Boundary points to the initial position
X1 = X1+repmat(mesh.X_CG,size(X1,1),1);

%Translating the center of small circles
mesh.X_CG_Small = mesh.X_CG_Small + mesh.X_CG;

%Creating the X2 Points of the Particle Elements
X2 = X1;
X2(1:end-1,:) = X1(2:end,:);
X2(end,:) = X1(1,:);



mesh.X1_P = X1;

mesh.X2_P = X2;

%Translating the real XCG
%NOTE THAT THE INITIAL XCG VALUE WILL NOW BE THE
%ACTUAL XCG, THEREFORE THE INITIAL XCG WILL NOT BE STORED ANYWHERE
% mesh.realX_CG = [x, y];
mesh.X_CG(1) = mesh.X_CG(1) + x ;
mesh.X_CG(2) = mesh.X_CG(2) + y ;

% %Putting the first Real Cent. Of Gravity point to the matrix which will
% %store all real center of gravities
% mesh.AllrealX_CG(1,:) = mesh.realX_CG;



%Making sure that this is the first and last time the points are created.
mesh.counter = mesh.counter + 1;

end

%At this point, if this was the first call of the function, the points are
%created and translated to their initial position. From now on, the points
%will only be translated and rotated according to the realXCG.

%Making sure that the translation and rotation in the iterations are not
%done in the initial creation of points




if mesh.iter
    
    %Part 1:                            -R O T A T I O N-
    
    %In order to rotate all points according to the realXCG, we first
    %translate the realXCG which was before the iteration to the origin, rotate, then
    %move the points to the current realXCG 
    v =  values.allX_CG(:,:,mesh.counter);
    x = v(1);
    y = v(2);
   
    %Moving the points to the origin
    mesh.X1_P(:,1) = mesh.X1_P(:,1) - x;
    mesh.X1_P(:,2) = mesh.X1_P(:,2) - y;
    
    
    %Rotating the X1 points
    [angle_i, dist_i] = cart2pol(mesh.X1_P(:,1) , mesh.X1_P(:,2) );
    
    %Adding the rotation amount (IN RADIANS)
    angle_i = angle_i + theta;
    
    %New rotated points
    [mesh.X1_P(:,1), mesh.X1_P(:,2)] = pol2cart(angle_i,dist_i);
    
    %Moving the points back to their places
    mesh.X1_P(:,1) = mesh.X1_P(:,1) + x;
    mesh.X1_P(:,2) = mesh.X1_P(:,2) + y;
   
    
    %Part 2:                        -T R A N S L A T I O N-
    %X1 points
    
    %translation vector
    v = mesh.X_CG - values.allX_CG(:,:,mesh.counter);
    
    mesh.X1_P(:,1) = mesh.X1_P(:,1) + v(1);
    mesh.X1_P(:,2) = mesh.X1_P(:,2) + v(2);
    
    %Creating the X2 Points of the Particle Elements
    mesh.X2_P = mesh.X1_P;
    mesh.X2_P(1:end-1,:) = mesh.X1_P(2:end,:);
    mesh.X2_P(end,:) = mesh.X1_P(1,:);
    
    
%Updating the mesh counter
mesh.counter = mesh.counter + 1;

    

    
    
    
    
    
end
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
    }
    %     Matrix containing the indices of the last elements of each circle
mesh.lastE_P = mesh.N_E_P;
% 
% % Matrix containing the indices of the first elements of each circle
  mesh.firstE_P = 1;
    
    
    
    
    
    
    Finding the differences between the coordinates of the midpoints 
    of the elements on the particle and the coordinates of the center
    of gravity of the paticle
     
    S : abbreviation for Segment
    
    Tg : abbreviation for Tangent
    
    STg : N_Ep x 2 matrix
    %}
    mesh.lastE_P = mesh.N_E_P;
    mesh.firstE_P = 1;
    % Initializing the relative position matrix
    mesh.S_Tg = zeros(sum(mesh.N_E_P), 2);
    
    %{
    Using the same index matrix used for filling up the coordinates 
    matrices
    REAL CENT OF GRAVITY SUBSTITUTED BELOW
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
    mesh.S_N = zeros(length(mesh.S_Tg),2);
    % x coordinates
    mesh.S_N(:, 1) = - mesh.S_Tg(:, 2) ./ mesh.S_L;
    
    % y coordinates
    mesh.S_N(:, 2) = mesh.S_Tg(:, 1) ./ mesh.S_L;
    
    
    
    




end