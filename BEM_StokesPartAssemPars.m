function [Tg, E_L, E_N, XM] = BEM_StokesPartAssemPars(X1, X2)
    % Assembly Parameters
    
    
    %%
    %{
    ***********************************************************************
        Relative Position Matrix
    ***********************************************************************
    
    Differences between the coordinates of the second points of the 
    elements and the coordinates of the first points of the elements
     
    Tg : abbreviation for Tangent which represents the aforementioned
    difference
    
    Tg : N x 2 matrix
    %}
    
	Tg = X2 - X1;
    
    
    %%
    %{
    ***********************************************************************
        Length of each element
    ***********************************************************************
    %}
    
    %{
    Length of each element found from the Pythagorean Theorem
    
    EL : abbreviation for ELength representing the aforementioned length
    
    EL : N x 1 vector
    %}
	E_L = hypot(Tg(:, 1), Tg(:, 2));
    
    
    %%
    %{
    ***********************************************************************
        Outward-pointing normal unit vector for each element
    ***********************************************************************
    
    EN : abbreviation for ENormal which represents the aforementioned
    normal unit vector
    %}
    
    %{
    x components of the unit vectors found from applying a dot product
    logic
    %}
    E_N(:, 1) = Tg(:, 2) ./ E_L;
    
    %{
    y components of the unit vectors found from applying a dot product
    logic
    %}
	E_N(:, 2) = - Tg(:, 1) ./ E_L;
    
    
    %%
    %{
    ***********************************************************************
        Midpoints Matrix
    ***********************************************************************
    Matrix containing the coordinates of the midpoints of the elements
    
    M : abbreviation for Midpoint
    
    XM : N x 2 matrix
    %}
    
    %{
    Matrix containing the coordinates of the midpoints of all the 
    elements
    %}
    XM = (X1 + X2) / 2; 
end