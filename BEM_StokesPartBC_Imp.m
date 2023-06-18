function [H, G] = BEM_StokesPartBC_Imp(N_E, BC, H, G)
    % Rearranging the entries of the matrices G and H according to the
    % boundary condition type
    
    
    % Iterating through all of the entries to separate the unknowns
    for i = 1 : 2 * N_E
        
        %{
        Changing the columns which correspond to the entries with a
        Dirichlet boundary condition
        %}
        if BC(i, 1) == 0
            
            % Introducing a dummy variable for the exchange 
            dummy = G(:, i);
            
            %{
            Changing the entries of the matrix G and adding a negative in
            front since the entry is being transferred to the other side
            of the equation
            %}
            G(:, i) = - H(:, i);
            
            %{
            Changing the entries of the matrix H in the same way as for the
            matrix G
            %}
            H(:,i) = - dummy;
        end
    end
end