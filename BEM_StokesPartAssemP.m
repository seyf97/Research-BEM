function  assembly = BEM_StokesPartAssemP(mesh, assembly)
    % Calculation of the system matrices for the Boundary Value Problem
    % given in a rectangular domain having a cylinder at its center which
    % is going be to solved with the Boundary Element Method
    
    
    %%
    %{
    ***********************************************************************
        Matrix that will be multiplied by the matrix containing the
        traction components at each element on the particle surface

    The product will give the augmented resultant vector at CG of the 
    particle.

        fB = F * tP

        F : (3 * Np) x (2 * N_E_P_T)
    ***********************************************************************
    %}

    % Initializing the matrix
    assembly.F = zeros(3 * mesh.N_P, 2 * mesh.N_E_P_O);

    %
    for i = 1 : mesh.N_P

            %
            assembly.F(3 * i - 2, 2 * mesh.firstE_P(i) - 1 : 2 : 2 * mesh.lastE_P(i))...
                = mesh.E_L_P(mesh.firstE_P(i) : mesh.lastE_P(i));

            %
            assembly.F(3 * i - 1, 2 * mesh.firstE_P(i) : 2 : 2 * mesh.lastE_P(i))...
                = mesh.E_L_P(mesh.firstE_P(i) : mesh.lastE_P(i));

            %
            assembly.F(3 * i,  2 * mesh.firstE_P(i) - 1 : 2 : 2 * mesh.lastE_P(i)) =...
                - mesh.S_Tg(mesh.firstE_P(i) : mesh.lastE_P(i), 2);

            %
            assembly.F(3 * i, 2 * mesh.firstE_P(i) : 2 : 2 * mesh.lastE_P(i)) =...
                mesh.S_Tg(mesh.firstE_P(i) : mesh.lastE_P(i), 1);
    end
    
    
    % In case particulate flow will be observed
    if ~mesh.findForce
        %%
        %{
        ***********************************************************************
            Matrix that will be multiplied by the rigid-body motion parameters
            of the particle

        The product will give the velocities at the computational nodes of
        the elements on the particle surface.

            uP = M * uB

            M : (2 * N_E_P_T) x (3 * Np)
        ***********************************************************************
        %}

        % Initializing the matrix
        assembly.M = zeros(2 * mesh.N_E_P_O, 3 * mesh.N_P);

        %
        for i = 1 : mesh.N_P

            %
            assembly.M(2 * mesh.firstE_P(i) - 1 : 2 : 2 * mesh.lastE_P(i), 3 * i - 2)...
                = 1;

            %
            assembly.M(2 * mesh.firstE_P(i) : 2 : 2 * mesh.lastE_P(i), 3 * i - 2)...
                = 0;  

            %
            assembly.M(2 * mesh.firstE_P(i) : 2 : 2 * mesh.lastE_P(i), 3 * i - 1)...
                = 1;

            %
            assembly.M(2 * mesh.firstE_P(i) - 1 : 2 : 2 * mesh.lastE_P(i), 3 * i - 1)...
                = 0; 
            %
            assembly.M(2 * mesh.firstE_P(i) - 1 : 2 : 2 * mesh.lastE_P(i), 3 * i)...
                = mesh.S_L(mesh.firstE_P(i) : mesh.lastE_P(i)) .* mesh.S_N(mesh.firstE_P(i) :...
                mesh.lastE_P(i), 1);

            %
            assembly.M(2 * mesh.firstE_P(i) : 2 : 2 * mesh.lastE_P(i), 3 * i)...
                = mesh.S_L(mesh.firstE_P(i) : mesh.lastE_P(i)) .* mesh.S_N(mesh.firstE_P(i) :...
                mesh.lastE_P(i), 2);
        end
    end
    
    
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
    N_E = mesh.N_E;
    
    %
    mu = mesh.mu;
    
    % Element length vector doubled
    E_L([1 : 2 : 2 * N_E 2 : 2 : 2 * N_E], 1) =...
        [mesh.E_L; mesh.E_L];
    
    % Matrix containing unit normal vectors to the elements doubled
    E_N([1 : 2 : 2 * N_E 2 : 2 : 2 * N_E], :) =...
        [mesh.E_N; mesh.E_N];
    
    % Midpoint matrix doubled
    XM([1 : 2 : 2 * N_E 2 : 2 : 2 * N_E], :) =...
        [mesh.XM; mesh.XM];
    
    % Matrix containing pieces of tangent vectors to the elements doubled
    Tg([1 : 2 : 2 * N_E 2 : 2 : 2 * N_E], :) =...
        [mesh.Tg; mesh.Tg];
    

HP0 = zeros(2*N_E_0);
GP0 = zeros(2*N_E_0);

   
    % The load point which will travel the rows of the H and G matrices
    parfor ld = 2*N_E_0 : 2*mesh.N_E
    
    % Initializing the rows to be filled in every iteration
    %singlerow1 is either row 1(1) or 2(1) 
    %and singlerow2 is either row 1(2) or 2(2) (column)
    HsingleRow1 = zeros(N_E_0,1);
    HsingleRow2 = HsingleRow1;
    
    GsingleRow1 = zeros(N_E_0,1);
    GsingleRow2 = GsingleRow1;
    
    %%
    % Field points will be evaluated below. Every field point will be
    % multiplied with the gauss points and weights in every iteration until
    % 12. The number 12 == number of abscissa points for the Gauss Quad. In
    % this assembly code, 12 points was observed to be sufficient for the
    % accuracy.
    
    % Note that until i=12, the values aren't completed yet! 
    for i=1:12
    
        % Summing all field points in the variables below
    x_i1 = XM(1:2*N_E_0,1) + Tg(1:2*N_E_0,1)/2*abscissa(i); %#ok<*PFBNS>
    x_i2 = XM(1:2*N_E_0,2) + Tg(1:2*N_E_0,2)/2*abscissa(i);
    
    rx = (x_i1 - XM(ld,1)).^2;
    ry = (x_i2 - XM(ld,2)).^2;
    
    r = sqrt( rx + ry );
    
    r_i1 = ( XM(ld,1) - x_i1 )./r;
    r_i2 = ( XM(ld,2) - x_i2 )./r;
    
    delrdeln = E_N(1:2*N_E_0,1) .* r_i1 + E_N(1:2*N_E_0,2) .* r_i2;
    
    % -----  Now we will find the mini matrix elements of the G and H matrices -----
    
    %     ---------- If the load point, or row number, is on the EVEN numbers  -------
    
            if mod(ld,2)==0
            
            %[22] Non Diagonal Elements
                   %H matrix
                  quadFuncH22 = r_i2(2:2:end) .* r_i2(2:2:end) .* delrdeln(2:2:end) ./ r(2:2:end);
    HsingleRow2 = HsingleRow2 + quadFuncH22 * weight(i);
            
                   %G matrix
                  quadFuncG22 = - log( r(2:2:end) ) + r_i2(2:2:end).^2;
    GsingleRow2 = GsingleRow2 + quadFuncG22 * weight(i);
            
            
            
       % ------------ If the load point, or row number, is on the ODD numbers -----------------
            elseif mod(ld,2)==1
            
            % [12] Non diagonal entries
            
            %H Matrix
            quadFuncH12 = r_i1(2:2:end) .* r_i2(2:2:end) .* delrdeln(2:2:end) ./ r(2:2:end);
     HsingleRow2 = HsingleRow2 + quadFuncH12 * weight(i);
     
            %G Matrix
            quadFuncG12 = r_i1(2:2:end) .* r_i2(2:2:end);
     GsingleRow2 = GsingleRow2 + quadFuncG12 * weight(i);

                   
            % [11] Non Diagonal Elements
            
            %H Matrix
            quadFuncH11 = r_i1(1:2:end) .* r_i1(1:2:end) .* delrdeln(1:2:end) ./ r(1:2:end);
      HsingleRow1 = HsingleRow1 + quadFuncH11 * weight(i);
      
            %G Matrix
            quadFuncG11 = -log( r(1:2:end) ) + r_i1(1:2:end).^2;
      GsingleRow1 = GsingleRow1 + quadFuncG11 * weight(i);
        
      
              %End of if statement
            end
   
    %End of calculating a single row, advancing down.           
    end
    
    %Allocating the even and odd numbered elements to their specified
    %locations
    rowH = zeros(2*N_E_0,1);
    rowG = rowH;
    rowH(1:2:end) = HsingleRow1;
    rowH(2:2:end) = HsingleRow2;
    
    rowG(1:2:end) = GsingleRow1;
    rowG(2:2:end) = GsingleRow2;
%     
%     % The diagonal elements of the G matrix are supposed to be calculated
%     % more precisely due to singularity
%      %NOTE THAT WE CHANGED ALL FD VARIABLES TO LD SINCE WE WILL ONLY USE
%      %THE BELOW FUNCTION HANDLES FOR THE DIAGONAL ELEMENTS OF G
    
    
    
    %Multiplying with the rest of the function variables. Note that we
    %leave the diagonal elements of the H matrix to the last part.
    rowG = rowG .* (E_L(1:2*N_E_0)/2);
    rowH = rowH .* (E_L(1:2*N_E_0)/2);
    
    rowG = 1/(4*pi*mu)*rowG';
    rowH = 1/(pi) * rowH';
    
 
    HP0(ld,:) = rowH';
    GP0(ld,:) = rowG';
    
    
    end
    
HP0(2:2:end,1:2:end) = HP0(1:2:end,2:2:end);
GP0(2:2:end,1:2:end) = GP0(1:2:end,2:2:end);

HP0 = HP0(2*N_E_0+1:end,:);
GP0 = GP0(2*N_E_0+1:end,:);





%%

%Matrix containing PP and 0P
HPP = zeros(2*mesh.N_E);
GPP = HPP;


   % The load point which will travel the rows of the H and G matrices
    parfor ld = 1 : 2*mesh.N_E

            isMatrixPP = false;
            if ld > 2*N_E_0
                isMatrixPP = true;
            end


        % Initializing the rows to be filled in every iteration
        %singlerow1 is either row 1(1) or 2(1) 
        %and singlerow2 is either row 1(2) or 2(2) (column)
        HsingleRow1 = zeros(mesh.N_E_P_O,1);
        HsingleRow2 = HsingleRow1;

        GsingleRow1 = zeros(mesh.N_E_P_O,1);
        GsingleRow2 = GsingleRow1;

        %%
        % Field points will be evaluated below. Every field point will be
        % multiplied with the gauss points and weights in every iteration until
        % 12. The number 12 == number of abscissa points for the Gauss Quad. In
        % this assembly code, 12 points was observed to be sufficient for the
        % accuracy.

        % Note that until i=12, the values aren't completed yet! 
        for i=1:12

            % Summing all field points in the variables below
        x_i1 = XM((2*N_E_0+1):end , 1) + Tg((2*N_E_0+1):end , 1)/2*abscissa(i);
        x_i2 = XM((2*N_E_0+1):end , 2) + Tg((2*N_E_0+1):end , 2)/2*abscissa(i);

        rx = (x_i1 - XM(ld,1)).^2;
        ry = (x_i2 - XM(ld,2)).^2;

        r = sqrt( rx + ry );

        r_i1 = ( XM(ld,1) - x_i1 )./r;
        r_i2 = ( XM(ld,2) - x_i2 )./r;

        delrdeln = E_N((2*N_E_0+1):end , 1) .* r_i1 + E_N((2*N_E_0+1):end , 2) .* r_i2;

        % -----  Now we will find the mini matrix elements of the G and H matrices -----

        %     ---------- If the load point, or row number, is on the EVEN numbers  -------

                if mod(ld,2)==0
                    %[22] Non Diagonal Elements
                    %H matrix
                    quadFuncH22 = r_i2(2:2:end) .* r_i2(2:2:end) .* delrdeln(2:2:end) ./ r(2:2:end);
                    HsingleRow2 = HsingleRow2 + quadFuncH22 * weight(i);
                    
                    %G matrix
                    quadFuncG22 = - log( r(2:2:end) ) + r_i2(2:2:end).^2;
                    GsingleRow2 = GsingleRow2 + quadFuncG22 * weight(i);
                    
                % ------------ If the load point, or row number, is on the ODD numbers -----------------
                elseif mod(ld,2)==1

                % [12] Non diagonal entries

                %H Matrix
                quadFuncH12 = r_i1(2:2:end) .* r_i2(2:2:end) .* delrdeln(2:2:end) ./ r(2:2:end);
                HsingleRow2 = HsingleRow2 + quadFuncH12 * weight(i);

                %G Matrix
                quadFuncG12 = r_i1(2:2:end) .* r_i2(2:2:end);
                GsingleRow2 = GsingleRow2 + quadFuncG12 * weight(i);


                % [11] Non Diagonal Elements

                %H Matrix
                quadFuncH11 = r_i1(1:2:end) .* r_i1(1:2:end) .* delrdeln(1:2:end) ./ r(1:2:end);
                HsingleRow1 = HsingleRow1 + quadFuncH11 * weight(i);

                %G Matrix
                quadFuncG11 = -log( r(1:2:end) ) + r_i1(1:2:end).^2;
                GsingleRow1 = GsingleRow1 + quadFuncG11 * weight(i);


                %End of if statement
                end

        end

        %Allocating the even and odd numbered elements to their specified
        %locations
        rowH = zeros(2*mesh.N_E_P_O,1);
        rowG = rowH;
        rowH(1:2:end) = HsingleRow1;
        rowH(2:2:end) = HsingleRow2;

        rowG(1:2:end) = GsingleRow1;
        rowG(2:2:end) = GsingleRow2;

        rowH = [zeros(2*N_E_0,1) ; rowH]; %#ok<*AGROW>
        rowG = [zeros(2*N_E_0,1) ; rowG];

        % The diagonal elements of the G matrix are supposed to be calculated
        % more precisely due to singularity
         %NOTE THAT WE CHANGED ALL FD VARIABLES TO LD SINCE WE WILL ONLY USE
         %THE BELOW FUNCTION HANDLES FOR THE DIAGONAL ELEMENTS OF G

                x_i = {@(t) XM(ld, 1) + Tg(ld, 1) / 2 .* t,...
                       @(t) XM(ld, 2) + Tg(ld, 2) / 2 .* t};
                %
                rxx = @(t) (x_i{1}(t) - XM(ld, 1)) .^ 2;
                ryy = @(t) (x_i{2}(t) - XM(ld, 2)) .^ 2;

                rr = @(t) sqrt(rxx(t) + ryy(t));
                %
                r_i = {@(t) (XM(ld, 1) - x_i{1}(t)) ./ rr(t),...
                       @(t) (XM(ld, 2) - x_i{2}(t)) ./ rr(t)};



        if isMatrixPP

            %Even diagonal elements
            if mod(ld,2) == 0
                quadFuncG = @(t) - log( rr(t) ) + r_i{2}(t) .^ 2;
                rowG(ld) = integral( quadFuncG , -1 , 1);

            %Odd diagonal elements
            else
                quadFuncG = @(t) - log( rr(t) ) + r_i{1}(t) .^ 2;
                rowG(ld) = integral( quadFuncG , -1 , 1);
            end

        end




        %Multiplying with the rest of the function variables. Note that we
        %leave the diagonal elements of the H matrix to the last part.
        rowG = rowG .* (E_L/2);
        rowH = rowH .* (E_L/2);

        rowG = 1/(4*pi*mu)*rowG';
        rowH = 1/(pi) * rowH';


        %Making sure that the diagonal elements of H are 0.5
        if isMatrixPP
            rowH(ld) = 0.5;
        end

        HPP(ld,:) = rowH';
        GPP(ld,:) = rowG';


    end

    HPP(2:2:end,1:2:end) = HPP(1:2:end,2:2:end);
    GPP(2:2:end,1:2:end) = GPP(1:2:end,2:2:end);

    H0P = HPP(1:2*N_E_0,(2*N_E_0+1):end);
    G0P = GPP(1:2*N_E_0,(2*N_E_0+1):end);

    HPP = HPP((2*N_E_0+1):end,(2*N_E_0+1):end);
    GPP = GPP((2*N_E_0+1):end,(2*N_E_0+1):end);

    
    %%
    %
    assembly.H = [assembly.H00 H0P;HP0 HPP];

    %
    assembly.G = [assembly.G00 G0P;GP0 GPP];
    
    %
    assembly.H0P = H0P;
    
    %
    assembly.HP0 = HP0;
    
    %
    assembly.HPP = HPP;
    
    %
    assembly.G0P = G0P;
    
    %
    assembly.GP0 = GP0;
    
    %
    assembly.GPP = GPP;
    
     
end