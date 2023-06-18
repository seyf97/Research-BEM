function BEM_StokesPartVel0Plot(mesh, values)
    % Plotting the velocity on the channel
    % 
    
    % Plotting the velocity on the channel
    quiver(mesh.XM_0(:, 1), mesh.XM_0(:, 2), values.u0(1 : 2 : end),...
        values.u0(2 : 2 : end))