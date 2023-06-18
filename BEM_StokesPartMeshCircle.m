function [X1, X2] =...
        BEM_StokesPartMeshCircle(Rad,X_CG,rotAng,mesh)
% Discretization of a circular boundary
% 
% Note that the angles are in radians.
    
    
    %%
    %{
    ***********************************************************************
        Coordinates of the first points of the elements on a circle
    
    X1 : N_Ep x 2 matrix
    ***********************************************************************
    %}
    
    
    
    
    % The function handle for calculating the x coordinates
    fx = @(argAng) Rad * cos(argAng + rotAng);
    
    % The function handle for calculating the y coordinates
    fy = @(argAng) Rad * sin(argAng + rotAng);
    
    % Full Angle: 360 degrees or 2 * pi
    fullAng = 2 * pi;
    
    % Sector angle pertaining to each element from the full circle
    secAng = fullAng / mesh.N_E_P;
    
    % Angle order in counter-clockwise direction
    angRange = 0 : secAng : fullAng;
    
    % Forming the points
    X1 = [fx(angRange)' fy(angRange)'];
    
    % Deleting the last point since it overlaps the 1st one
    X1(end, :) = [];
    
    % Making the order clockwise due to the convention in use
    X1(2 : end, :) = X1(end : - 1 : 2, :);
    
    % Bringing the x coordinates to the center of the rectangle
    X1(:, 1) = X1(:, 1) + X_CG(1);
    
    % Bringing the y coordinates to the center of the rectangle
    X1(:, 2) = X1(:, 2) + X_CG(2);
    
    
    %%
    %{
    ***********************************************************************
        Coordinates of the second points of the elements of the circle
    
    X2 : abbreviation for 2nd coordinates
    
    X2 : N_Ep x 2 matrix
    ***********************************************************************
    %}
    
    %{
    Preparing the elements of the circular boundary to be added to the
    second points vector
    
    X2 : N_Ep x 2 matrix
    %}
    X2 = X1;
    
    %
    X2 = [X2; X1(1, :)];
    
    %
    X2(1, :) = [];
    
    
if ~mesh.iter
%NUMBER OF ELEMENTS
mesh.N_E_P = mesh.N_E_P_O;

%Matrix containing the indices of the last elements of each circle
mesh.lastE_P = mesh.N_E_P;

% Matrix containing the indices of the first elements of each circle
mesh.firstE_P = mesh.lastE_P - mesh.N_E_P + 1;
end



end