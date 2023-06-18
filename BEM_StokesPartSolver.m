function [u, t] = BEM_StokesPartSolver(N_E, BC, H, G)
    %This function stores u and t values in seperate vectors.
    
    %
    u = zeros(2 * N_E, 1);
    
    %
    t = zeros(2 * N_E, 1);
    
    %{
    Solves the equation of the form A*x=b with the matrix A being the
    rearranged H matrix and the vector b being the matrix multiplication
    of the G matrix with the vector containing the boundary values
    %}

    %{
    The right-hand side vector of Ax=b as defined in the description above
    %}
    b = G * BC(:, 2);
    
    % The solution, namely the unknown boundary values
    x = H \ b;
    
    %{
    This loop is for checking the boundary condition and generate u and q
    vectors
    %}
    for i = 1 : 2 * N_E
        
        %{
        If the boundary condition is drichlet, then u vector stores the 
        actual value
        %}
        if BC(i, 1) == 0
            %
            u(i) = BC(i, 2);
            
            %
            t(i) = x(i);
            
        %If it is not, then u vector contains the value comes from x vector
        else
            %
            u(i) = x(i);
            
            %
            t(i) = BC(i, 2);
        end
    end
    
end