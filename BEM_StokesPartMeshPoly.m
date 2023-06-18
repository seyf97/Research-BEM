function [X1_0, X2_0] = BEM_StokesPartMeshPoly...
    (N_S_Poly, N_E_Poly, N_E_0, firstE_0, lastE_0, X_Vc)
    % Forming the first and second points of the elements at the side of
    % the polygon
    
    %{
    Adding the coordinates of the first vertex to the end of the matrix
    to form a closed polygon which is necessary for the rest of the
    calculation below.
    %}
    X_Vc = [X_Vc; X_Vc(1, :)];
    
    % Initializing the matrix containing the 1st point coordinates
    X1_0 = zeros(N_E_0, 2);
    
    % Initializing the matrix containing the 2nd point coordinates
    X2_0 = zeros(N_E_0, 2);
    
    %{
    Forming the first and second points of the elements at a single side of
    the polygon
    %}
    for i = 1 : N_S_Poly
        [X1_0(firstE_0(i) : lastE_0(i), :),...
         X2_0(firstE_0(i) : lastE_0(i), :)] = BEM_StokesPartMeshLine...
            (X_Vc(i, :), X_Vc(i + 1, :), N_E_Poly(i));
    end
end