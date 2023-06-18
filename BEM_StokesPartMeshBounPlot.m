function BEM_StokesPartMeshBounPlot(mesh)
    % Plotting the elements on rectangular boundary
    
    % In case the first update has been made
    if ishandle(mesh.figNumBounMesh)
        %{
        Removing the figure with the specified number to avoid overwriting
        the same figure
        %}
        close(mesh.figNumBounMesh)
    end
    
    % Prompting the Boundary Figure
    figure(mesh.figNumBounMesh)
    
    % Channel Plot
    BEM_StokesPartMeshBounPlot0(mesh)
    
    % Plotting the particles
    plot(mesh.XM_P([1 : end 1], 1), mesh.XM_P([1 : end 1], 2),...
        'r', 'LineWidth', 1)
    
    % Marking the first node of each particle for illustrating the rotation
    scatter(mesh.XM_P(mesh.firstE_P, 1),mesh.XM_P(mesh.firstE_P, 2),...
        'b' , 'filled', 'LineWidth', 0.5)
    
    % In case the channel is rectangular
    if mesh.isSRC
        % Adjusting the aspect ratio
        pbaspect([mesh.Lx / mesh.Ly 1 1])
    end

    % Making the figure full screen programmatically
    set(gcf, 'Position', get(0, 'Screensize'));
    
    
    % Showing the particle index
    for i = 1 : mesh.N_P
        text(mesh.X_CG(i, 1), mesh.X_CG(i, 2),...
          ['   ' num2str(i)], 'HorizontalAlignment', 'left', 'FontSize', 7)
    end
end