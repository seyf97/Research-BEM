function BEM_StokesPartVelDomPlot(mesh, values)
    % Velocity Field Plot
    %
    
    % Checking if the figure exists
    if ishandle(mesh.figNumVelField)
        %{
        Removing the figure with the specified number to avoid overwriting
        the same figure
        %}
        close(mesh.figNumVelField)
    end
    
    % Opening the appropriate figure
    figure(mesh.figNumVelField)
    
    % Holding the figure
    hold on
    
    
    % Drawing lines between the vertices of the channel/polygon
    plot(mesh.X_Vc([1 : end 1], 1), mesh.X_Vc([1 : end 1], 2),...
        'LineWidth', 2.5)
    
    % Plotting the particles
    plot(mesh.XM_P([1 : end 1], 1), mesh.XM_P([1 : end 1], 2),...
        'r', 'LineWidth', 1)
    
    % Marking the first node of each particle for illustrating the rotation
    scatter(mesh.XM_P(mesh.firstE_P, 1),mesh.XM_P(mesh.firstE_P, 2),...
        'b' , 'filled', 'LineWidth', 0.5)
    
    
    quiver(mesh.X_D(:, 1), mesh.X_D(:, 2), values.uD(1 : 2 : end),...
                                 values.uD(2 : 2 : end), 'k')
    
    %
    if mesh.isSRC
        %
        pbaspect([mesh.Lx / mesh.Ly 1 1])
    end

    % Making the figure full screen programmatically
    set(gcf, 'Position', get(0, 'Screensize'));
                             
    % Releasing the figure that was held in function VelDomPlot
    hold off
end