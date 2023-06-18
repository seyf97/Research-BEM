function BEM_StokesPartMeshBounPlot0(mesh)
    % Drawing the channel
    
    
    % Drawing lines between the vertices of the channel/polygon
    plot(mesh.X_Vc([1 : end 1], 1), mesh.X_Vc([1 : end 1], 2),...
        'LineWidth', 2.5)
    
    %
    grid on
    
    %
    hold on
    
    %
    if mesh.isSRC
        %
        pbaspect([mesh.Lx / mesh.Ly 1 1])
    end
    
    % Making the current figure full screen programmatically
    set(gcf, 'Position', get(0, 'Screensize'));
end