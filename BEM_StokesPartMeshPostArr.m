function [X1_Post, X2_Post, X_CG_Post, N_Post, N_E_Post_O] =...
        BEM_StokesPartMeshPostArr...
        (offsetx, offsety, postRad, shifty,...
        postArrLx, postArrLy, X_CG_1, N_E_Post)
    % Mesh Post Array

    %{
    500 is max number of posts when slipy is zero. Thus, we initialize
    the matrix storing the first points with the maximum possible number
    of elements and afterwards, dlete the trailing zeros if necessary
    %}
    X1_Post = zeros(500 * N_E_Post, 2);
    
    %{
    Initializing the matrix storing the center coordinates of the post
    array
    
    We will delete the trailing zeros if necessary after filing up the 
    matrix.
    %}
    X_CG_Post = zeros(500, 2);

    %
    inity = X_CG_1(2);

    %
    X_CG = X_CG_1;
    
    %
    X_CG_Post(1, :) = X_CG;

    %
    i = 1;

    %%
	%{
    ***********************************************************************
        1st Point Coordinates
    ***********************************************************************
    
    We calculate the coordinates of the first points of the elements.
    %}
    while 1
        
        %
        if X_CG(1) + postRad >= postArrLx
            break

        %
        elseif X_CG(2) + postRad >= postArrLy
            
            %
            if X_CG_1(2) + shifty >= inity + offsety

                %
                X_CG_1(1) = X_CG_1(1) + offsetx;

                %
                X_CG_1(2) = inity + shifty;

                %
                X_CG = X_CG_1;

            else
                %
                X_CG_1(1) = X_CG_1(1) + offsetx;

                %
                X_CG_1(2) = X_CG_1(2) + shifty;

                %
                X_CG = X_CG_1;
            end

        %
        else
            %
            X1_Post((i - 1) * N_E_Post + 1 : i * N_E_Post, :) =...
                BEM_StokesPartMeshCircle...
                (postRad, N_E_Post, X_CG, 0);

            %
            i = i + 1;

            %
            X_CG(2) = X_CG(2) + offsety;

        end
            
        %
        X_CG_Post(i, :) = X_CG;
    end

    % Deleting trailing zeros
    X1_Post = X1_Post(1 : find(X1_Post(:, 1), 1, 'last'), :);

    % Deleting trailing zeros
    X_CG_Post = X_CG_Post(1 : find(X_CG_Post(:, 1), 1, 'last'), :);
    
    %{
    Deleting the last element of center points because after all the
    posts are done i is still incremented and X_CG gets the last element
    out of the channel.
    %}
    X_CG_Post(end, :) = [];
    
    % Number of posts
    N_Post = size(X1_Post, 1) / N_E_Post;
    
    % Total Number of Elements in the post array
    N_E_Post_O = N_Post * N_E_Post;
    
    
    %%
	%{
    ***********************************************************************
        2nd Point Coordinates
    ***********************************************************************
    
    We calculate the coordinates of the second points of the elements.
    %}

    %
    X2_Post = [X1_Post; [0 0]];
    
    %
    for i = N_Post : - 1 : 1
        
        %
        X2_Post(i * N_E_Post + 1, :) = X2_Post((i - 1) * N_E_Post + 1, :);
    end
    
    %
    X2_Post(1, :) = [];
end