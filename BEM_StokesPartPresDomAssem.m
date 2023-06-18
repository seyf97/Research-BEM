function assembly = BEM_StokesPartPresDomAssem(mesh, assembly)
    % Calculation of the system matrices for finding the pressure field in
    % the domain
    % 
    % 
    % Note that the field points are the nodes of the elements at the
    % boundaries for post-processing.
    % 
    % 
    % mu : Dynamic Viscosity of the fluid [Pa.s]
    % 
    % N_E : Number of elements at the boundaries [1 x 1]
    % 
    % N_D : Number of points in the domain [1 x 1]
    % 
    % X_D : Matrix containing the positions of the nodes in the domain 
    %       [N_D x 2]
    % 
    % XM : Matrix containing the positions of the nodes at the boundary 
    %       [N_E x 2]
    % 
    % E_L : Vector containing the lengths of the elements at the boundary
    %       [N_E x 1]
    % 
    % E_N : Matrix containing the components of the outward unit normal 
    %       vectors of the elements at the boundary [N_E x 2]
    % 
    % Tg : Matrix containing the relative positions of the elements at the
    %      boundary [N_E x 2]
    
    
    %%
    %{
    ***********************************************************************
    Gauss Points
    ***********************************************************************
    
    Here, we use 12 Gauss points for integral evaluation which has so far 
    sufficed for the precision demand of the program.
    %}
    
    %{
    Weights
    
    These are multiplied by the integrands.
    %}
    weight = [0.2491470458134028 0.2491470458134028 0.2334925365383548...
              0.2334925365383548 0.2031674267230659 0.2031674267230659...
              0.1600783285433462 0.1600783285433462 0.1069393259953184...
              0.1069393259953184 0.0471753363865118 0.0471753363865118];
    
    %{
    Abscissas
              
    These are plugged into the integrand functions. As pointed out above, 
    the results are multiplied by the weights and summed up to calculate
    an approximate value for the integral at hand.
    %}
    abscissa =...
         [-0.1252334085114689  0.1252334085114689 -0.3678314989981802...
           0.3678314989981802 -0.5873179542866175  0.5873179542866175...
          -0.7699026741943047  0.7699026741943047 -0.9041172563704749...
           0.9041172563704749 -0.9815606342467192  0.9815606342467192];

    
    %%
    %{
    ***********************************************************************
        System Matrix Assembly
    ***********************************************************************
    %}
    
    %{
    Initializing the matrix containing the fundamental solutions of
    pressure across the elements.
    %}
    assembly.Q = zeros(mesh.N_D, 2 * mesh.N_E);
    
    %{
    Initializing the matrix containing the fundamental solutions of
    traction across the elements
    %}
    assembly.X = assembly.Q;
    
    %{
    Forming the load-point for loop only travelling the load points in
    the rectangular portion
    %}
    for ld = 1 : mesh.N_D
        %{
            Initializing the rows to be filled in every iteration

        singlerow : row 1(1)
 
        singlerow : row 1(2) 

        () indicates column number.
        %}

        %
        QsingleRow1 = zeros(mesh.N_E, 1);
        
        %
        XsingleRow1 = QsingleRow1;

        %
        QsingleRow2 = QsingleRow1;
        
        %
        XsingleRow2 = QsingleRow2;
        
        %{
        Forming the field-point for loop only travelling the field
        points in the rectangular portion
        %}
        for i = 1 : length(abscissa)
            % Summing all field points in the variables below
            x_i = {mesh.XM(:, 1) + mesh.Tg(:, 1) / 2 * abscissa(i),...
                   mesh.XM(:, 2) + mesh.Tg(:, 2) / 2 * abscissa(i)};
            
            r = hypot(x_i{1} - mesh.X_D(ld, 1), x_i{2} - mesh.X_D(ld, 2));
            
            r_i = {(mesh.X_D(ld, 1) - x_i{1}) ./ r,...
                   (mesh.X_D(ld, 2) - x_i{2}) ./ r};
            
            r_n = mesh.E_N(:, 1) .* r_i{1} + mesh.E_N(:, 2) .* r_i{2};
            
            
            %%
            %{
            ***************************************************************
            Q matrix
            ***************************************************************
            %}
            
            %[11]
            fundSolQ1 = r_i{1} ./ r;
            
            %
            QsingleRow1 = QsingleRow1 + fundSolQ1 * weight(i);
            
            %[12]
            fundSolQ2 = r_i{2} ./ r;
            
            %
            QsingleRow2 = QsingleRow2 + fundSolQ2 * weight(i);
            
            
            %%
            %{
            ***************************************************************
            X matrix
            ***************************************************************
            %}
            
            %[11]
            fundSolX1 = (mesh.E_N(:, 1) - 2 * r_i{1} .* r_n) ./ r .^ 2;
            
            %
            XsingleRow1 = XsingleRow1 + fundSolX1 * weight(i);
            
            %[12]
            fundSolX2 = (mesh.E_N(:, 2) - 2 * r_i{2} .* r_n) ./ r .^ 2;
            
            %
            XsingleRow2 = XsingleRow2 + fundSolX2 * weight(i);
            
        %{
        We are done with finding the entries for the load point ld once we 
        are out of this inner loop
        %}
        end
        
        %Allocating the even and odd numbered elements to their specified
        %locations
        rowQ = zeros(2 * mesh.N_E, 1);
        rowX = rowQ;
        
        rowQ(1:2:end) = QsingleRow1 .* (mesh.E_L / 2);
        rowQ(2:2:end) = QsingleRow2 .* (mesh.E_L / 2);

        rowX(1:2:end) = XsingleRow1 .* (mesh.E_L / 2);
        rowX(2:2:end) = XsingleRow2 .* (mesh.E_L / 2);

        
        %Multiplying with the rest of the function variables.
        assembly.Q(ld, :) = 1 / (2 * pi) * rowQ;
        assembly.X(ld, :) = mesh.mu / pi * rowX;
        
    %{
    We are done with finding the entries for all the load points once we 
    are out of this outer loop
    %}
    end
end