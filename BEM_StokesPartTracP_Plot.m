function BEM_StokesPartTracP_Plot(mesh, values)
    % Plotting the traction on the particle
    % 
    % Also, negating the tractions since they are positive when the stress 
    % is directed inside the system. However, since we desire to see the 
    % effect of the flow on the particle we are reversing the direction for 
    % illustration.
    
    
    quiver(mesh.XM_P(:, 1), mesh.XM_P(:, 2), - values.tP(1 : 2 : end),...
                                             - values.tP(2 : 2 : end))