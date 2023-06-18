function assembly = BEM_StokesPartVelDomAssem(mesh, assembly)
    % Calculation of the system matrices for the Boundary Value Problem
    % given in a rectangular domain having a cylinder at its center which
    % is going be to solved with the Boundary Element Method
    
    
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
    Preparation of Assembly Parameters
    ***********************************************************************
    
    Doubling the size of the assembly parameters in order to utilize them 
    for the calculations of the system matrices below.
    %}
    
    % Element length matrix doubled
    E_L([1 : 2 : 2 * mesh.N_E 2 : 2 : 2 * mesh.N_E], 1) = [mesh.E_L; mesh.E_L];
    
    % Matrix containing unit normal vectors to the elements doubled
    E_N([1 : 2 : 2 * mesh.N_E 2 : 2 : 2 * mesh.N_E], :) = [mesh.E_N; mesh.E_N];
   
    % Midpoint matrix doubled
    XM([1 : 2 : 2 * mesh.N_E 2 : 2 : 2 * mesh.N_E], :) = [mesh.XM; mesh.XM];
   
    % Load Point matrix doubled
    X_D([1 : 2 : 2 * mesh.N_D 2 : 2 : 2 * mesh.N_D], :) = [mesh.X_D; mesh.X_D];
    
    % Matrix containing pieces of tangent vectors to the elements doubled
    Tg([1 : 2 : 2 * mesh.N_E 2 : 2 : 2 * mesh.N_E], :) = [mesh.Tg; mesh.Tg];
    
    
    %%
    %{
    ***********************************************************************
        System Matrix Assembly
    ***********************************************************************
    %}
    
    % Initializing the H00 matrix
    assembly.P = zeros(2 * mesh.N_D, 2 * mesh.N_E);
    
    % Initializing the G00 matrix
    assembly.W = assembly.P;
    
    
    %%
    % Load points which will travel the rows of the H and G matrices
    for ld = 1 : 2 * mesh.N_D
        %{
            Initializing the rows to be filled in every iteration

        singlerow1 : either row 1(1) or 2(1) 

        singlerow2 : either row 1(2) or 2(2) 

        () indicates column number.
        %}

        %
        PsingleRow1 = zeros(mesh.N_E, 1);
        
        %
        PsingleRow2 = PsingleRow1;
        
        %
        WsingleRow1 = PsingleRow1;
        
        %
        WsingleRow2 = WsingleRow1;
        
        
        %%
        % Field points will be evaluated below. Every field point will be
        % multiplied with the gauss points and weights in every iteration until
        % 12. The number 12 == number of abscissa points for the Gauss Quad.
        
        % Note that until i = 12, the values aren't completed yet! 
        for i = 1 : 12
            % Summing all field points in the variables below
            x_i1 = XM(:, 1) + Tg(:, 1) / 2 * abscissa(i);
            
            x_i2 = XM(:, 2) + Tg(:, 2) / 2 * abscissa(i);
            
            rx = (x_i1 - X_D(ld, 1)) .^ 2;
            ry = (x_i2 - X_D(ld, 2)) .^ 2;
            
            r = sqrt(rx + ry);
            
            r_i1 = (X_D(ld, 1) - x_i1) ./ r;
            r_i2 = (X_D(ld, 2) - x_i2) ./ r;
            
            r_n = E_N(:,1) .* r_i1 + E_N(:,2) .* r_i2;
            
            % Now we will find the mini matrix elements of the G and H matrices
            
            % If the load point, or row number, is on the EVEN numbers
            
            if mod(ld,2)==0
                %[22] Non Diagonal Elements
                %H matrix
                quadFuncP22 = r_i2(2:2:end) .* r_i2(2:2:end) .* r_n(2:2:end) ./ r(2:2:end);
                PsingleRow2 = PsingleRow2 + quadFuncP22 * weight(i);
                
                %W Matrix
                quadFuncW22 = - log( r(2:2:end) ) + r_i2(2:2:end).^2;
                WsingleRow2 = WsingleRow2 + quadFuncW22 * weight(i);
            % If the load point, or row number, is on the ODD numbers
            elseif mod(ld,2)==1
                % [12] Non diagonal entries
                %P Matrix
                quadFuncP12 = r_i1(2:2:end) .* r_i2(2:2:end) .* r_n(2:2:end) ./ r(2:2:end);
                PsingleRow2 = PsingleRow2 + quadFuncP12 * weight(i);

                %W Matrix
                quadFuncW12 = r_i1(2:2:end) .* r_i2(2:2:end);
                WsingleRow2 = WsingleRow2 + quadFuncW12 * weight(i);


                % [11] Non Diagonal Elements
                %P Matrix
                quadFuncP11 = r_i1(1:2:end) .* r_i1(1:2:end) .* r_n(1:2:end) ./ r(1:2:end);
                PsingleRow1 = PsingleRow1 + quadFuncP11 * weight(i);

                %W Matrix
                quadFuncW11 = -log( r(1:2:end) ) + r_i1(1:2:end).^2;
                WsingleRow1 = WsingleRow1 + quadFuncW11 * weight(i);
            end
        end

        %Allocating the even and odd numbered elements to their specified
        %locations
        rowP = zeros(2*mesh.N_E,1);
        rowW = rowP;
        
        rowP(1:2:end) = PsingleRow1;
        rowP(2:2:end) = PsingleRow2;

        rowW(1:2:end) = WsingleRow1;
        rowW(2:2:end) = WsingleRow2;

        
        %Multiplying with the rest of the function variables. 
        rowW = rowW .* (E_L/2);
        rowP = rowP .* (E_L/2);

        rowW = 1/(4*pi*mesh.mu)*rowW';
        rowP = 1/(pi) * rowP';


        assembly.P(ld, :) = rowP';
        assembly.W(ld, :) = rowW';
    end
    
    assembly.P(2:2:end,1:2:end) = assembly.P(1:2:end,2:2:end);
    assembly.W(2:2:end,1:2:end) = assembly.W(1:2:end,2:2:end);
end