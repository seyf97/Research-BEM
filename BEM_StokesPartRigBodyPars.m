function values = BEM_StokesPartRigBodyPars(assembly, values)
    %{
    Solves the equation of the form A*x=b with the matrix A being the
    rearranged H matrix and the vector b being the matrix multiplication
    of the G matrix with the vector containing the boundary values

        uB : 3 x 1
    %}

    %{
    The right-hand side vector of A * x = b as defined in the description
    above
    %}
    b = - assembly.R * values.u0;
    
    % The solution, namely the rigid-body motion parameters
    values.uB = assembly.T \ b;