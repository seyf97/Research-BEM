function assembly = BEM_StokesPartAssem0(mesh)
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
    
    %
    N_E_0 = mesh.N_E_0;
    
    %
    mu = mesh.mu;
    
    % Element length vector doubled
    E_L_0([1 : 2 : 2 * mesh.N_E_0 2 : 2 : 2 * mesh.N_E_0], 1) =...
        [mesh.E_L_0; mesh.E_L_0];
    
    % Matrix containing unit normal vectors to the elements doubled
    E_N_0([1 : 2 : 2 * mesh.N_E_0 2 : 2 : 2 * mesh.N_E_0], :) =...
        [mesh.E_N_0; mesh.E_N_0];
    
    % x components of E_N_0
    E_N_01 = E_N_0(:, 1);
    
    % y components of E_N_0
    E_N_02 = E_N_0(:, 2);
    
    % Midpoint matrix doubled
    XM_0([1 : 2 : 2 * mesh.N_E_0 2 : 2 : 2 * mesh.N_E_0], :) =...
        [mesh.XM_0; mesh.XM_0];
    
    % x components of XM_0
    XM_01 = XM_0(:, 1);
    
    % y components of XM_0
    XM_02 = XM_0(:, 2);
    
    % Matrix containing pieces of tangent vectors to the elements doubled
    Tg_0([1 : 2 : 2 * mesh.N_E_0 2 : 2 : 2 * mesh.N_E_0], :) =...
        [mesh.Tg_0; mesh.Tg_0];
    
    % x components of Tg_0
    Tg_01 = Tg_0(:, 1);
    
    % y components of Tg_0
    Tg_02 = Tg_0(:, 2);
    
    
    %%
    %{
    ***********************************************************************
        System Matrix Assembly
    ***********************************************************************
    %}
    
    % Initializing the H00 matrix
    H00 = zeros(2 * N_E_0);
    
    % Initializing the G00 matrix
    G00 = H00;
    
    
    %%
    % Load points which will travel the rows of the H and G matrices
    parfor ld = 1 : 2 * N_E_0
        %{
            Initializing the rows to be filled in every iteration

        singlerow1 : either row 1(1) or 2(1) 

        singlerow2 : either row 1(2) or 2(2) 

        () indicates column number.
        %}

        %
        HsingleRow1 = zeros(N_E_0, 1);
        
        %
        HsingleRow2 = HsingleRow1;
        
        %
        GsingleRow1 = HsingleRow1;
        
        %
        GsingleRow2 = GsingleRow1;
        
        
        %%
        % Field points will be evaluated below. Every field point will be
        % multiplied with the gauss points and weights in every iteration until
        % 12. The number 12 == number of abscissa points for the Gauss Quad.
        
        % Note that until i = 12, the values aren't completed yet! 
        for i = 1 : 12
            % Summing all field points in the variables below
            x_i = {XM_01 + Tg_01 / 2 * abscissa(i),...
                   XM_02 + Tg_02 / 2 * abscissa(i)}; %#ok<*PFBNS>
            
            r = hypot(XM_01(ld) - x_i{1}, XM_02(ld) - x_i{2});
            
            r_i = {(XM_01(ld) - x_i{1}) ./ r,...
            	   (XM_02(ld) - x_i{2}) ./ r};
            
            r_n = E_N_01 .* r_i{1} + E_N_02 .* r_i{2};
            
            % Now we will find the mini matrix elements of the G and H matrices
            
            % If the load point, or row number, is on the EVEN numbers
            
            if mod(ld, 2) == 0
                %[22] Non Diagonal Elements
                %H matrix
                fundSolH22 = r_i{2}(2 : 2 : end) .* r_i{2}(2 : 2 : end)...
                              .* r_n(2 : 2 : end) ./ r(2 : 2 : end);
                          
                HsingleRow2 = HsingleRow2 + fundSolH22 * weight(i);
                
                %G matrix
                fundSolG22 = - log(r(2 : 2 : end)) +...
                    r_i{2}(2 : 2 : end) .^ 2;
                
                GsingleRow2 = GsingleRow2 + fundSolG22 * weight(i);
            % If the load point, or row number, is on the ODD numbers
            elseif mod(ld, 2) == 1
                % [12] Non diagonal entries
                %H Matrix
                fundSolH12 = r_i{1}(2 : 2 : end) .* r_i{2}(2 : 2 : end)...
                             .* r_n(2 : 2 : end) ./ r(2 : 2 : end);
                         
                HsingleRow2 = HsingleRow2 + fundSolH12 * weight(i);

                %G Matrix
                fundSolG12 = r_i{1}(2 : 2 : end) .* r_i{2}(2 : 2 : end);
                
                GsingleRow2 = GsingleRow2 + fundSolG12 * weight(i);


                % [11] Non Diagonal Elements
                %H Matrix
                fundSolH11 = r_i{1}(1 : 2 : end) .* r_i{1}(1 : 2 : end)...
                             .* r_n(1 : 2 : end) ./ r(1 : 2 : end);
                         
                HsingleRow1 = HsingleRow1 + fundSolH11 * weight(i);

                %G Matrix
                fundSolG11 = -log(r(1 : 2 : end)) +...
                    r_i{1}(1 : 2 : end) .^ 2;
                
                GsingleRow1 = GsingleRow1 + fundSolG11 * weight(i);
            end
        end

        %Allocating the even and odd numbered elements to their specified
        %locations
        rowH = zeros(2 * N_E_0, 1);
        rowG = rowH;
        rowH(1 : 2 : end) = HsingleRow1;
        rowH(2 : 2 : end) = HsingleRow2;

        rowG(1 : 2 : end) = GsingleRow1;
        rowG(2 : 2 : end) = GsingleRow2;

        % The diagonal elements of the G matrix are supposed to be calculated
        % more precisely due to singularity
        %NOTE THAT WE CHANGED ALL FD VARIABLES TO LD SINCE WE WILL ONLY USE
        %THE BELOW FUNCTION HANDLES FOR THE DIAGONAL ELEMENTS OF G
        
        x_i = {@(t) XM_01(ld) + Tg_01(ld) / 2 .* t,...
               @(t) XM_02(ld) + Tg_02(ld) / 2 .* t};
        
        r = @(t) hypot(XM_01(ld) - x_i{1}(t), XM_02(ld) - x_i{2}(t));
        %
        r_i = {@(t) (XM_01(ld) - x_i{1}(t)) ./ r(t),...
               @(t) (XM_02(ld) - x_i{2}(t)) ./ r(t)};
        
        %Even diagonal elements
        if mod(ld,2) == 0
            fundSolG = @(t) - log(r(t)) + r_i{2}(t) .^ 2;
            
            rowG(ld) = integral(fundSolG , -1 , 1);

        %Odd diagonal elements
        else
            fundSolG = @(t) - log(r(t)) + r_i{1}(t) .^ 2;
            
            rowG(ld) = integral(fundSolG , -1 , 1);
        end

        %Multiplying with the rest of the function variables. Note that we
        %leave the diagonal elements of the H matrix to the last part.
        rowG = rowG .* (E_L_0 / 2);
        rowH = rowH .* (E_L_0 / 2);

        rowG = 1 / (4 * pi * mu) * rowG';
        rowH = 1 / pi * rowH';

        %Making sure that the diagonal elements of H are 0.5
        rowH(ld) = 0.5;


        H00(ld, :) = rowH';
        G00(ld, :) = rowG';
    end
    
    H00(2 : 2 : end, 1 : 2 : end) = H00(1 : 2 : end, 2 : 2 : end);
    G00(2 : 2 : end, 1 : 2 : end) = G00(1 : 2 : end, 2 : 2 : end);
    
    
    %%
    %{
    ***********************************************************************
        Inversing G00 matrix once here at the begining to save from the
        computational cost since it will be used throughout.
    ***********************************************************************
    %}
    assembly.G00Inv = eye(2 * mesh.N_E_0) / G00;
    
    
    %%
    %{
        Putting the variables inside the appropriate structures
    %}
    
    %
    assembly.H00 = H00;
    
    %
    assembly.G00 = G00;
    
    %
    
end