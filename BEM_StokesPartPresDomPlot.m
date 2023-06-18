function mesh = BEM_StokesPartPresDomPlot(mesh, assembly, values)
    % Demonstrating the pressure distribution inside the domain 
    %
    % We are plotting over the velocity field that was called before this
    % function in the function Drag
    
    % Checking if the figure exists
    if ishandle(mesh.figNumPresField)
        %{
        Removing the figure with the specified number to avoid overwriting
        the same figure
        %}
        close(mesh.figNumPresField)
    end
    
    % Opening the appropriate figure
    figure(mesh.figNumPresField)
    
    % Holding the figure
    hold on
    
    % Drawing lines between the vertices of the channel/polygon
    plot(mesh.X_Vc([1 : end 1], 1), mesh.X_Vc([1 : end 1], 2),...
        'LineWidth', 2.5)
    
    % Pressure Field visualized
    if mesh.isSRC && ~mesh.isDLD
        %
        x = mesh.XM_0(1 : mesh.N_E_0x, 1);
        
        %
        y = mesh.XM_0(end - mesh.N_E_0y + 1 : end, 2);
        
        %
        x = x(~(x > mesh.X_CG(1) - mesh.Rad & x < mesh.X_CG(1) + mesh.Rad));
        
        %
        y = y(~(y > mesh.X_CG(2) - mesh.Rad & y < mesh.X_CG(2) + mesh.Rad));
        mesh.x=x;mesh.y=y;
        %
        [xx, yy] = meshgrid(x, y);
        mesh.xx=xx;mesh.yy=yy;
        %
        mesh.X_D = [xx(:) yy(:)];
        
        %
        mesh.N_D = length(mesh.X_D);
        
        %
        [~, values] = BEM_StokesPartPresDom(mesh, assembly, values);
        
        %
        contourf(xx, yy, reshape(values.pD, length(y), length(x)))
    
        %
        if mesh.isSRC
            %
            pbaspect([mesh.Lx / mesh.Ly 1 1])
        end
        
        %
%         scatter3(xx(:), yy(:), values.pD)
    else
        %
        scatter3(mesh.X_D(:, 1), mesh.X_D(:, 2), values.pD);
        
        % Adjusting the viewpoint
        view([12 22])
    end
    
    % Interpolating in between the nodes
    shading interp
    
    % Showing the correspoding value range of the colors
    colorbar
    
    % Plotting the particles
    plot(mesh.XM_P([1 : end 1], 1), mesh.XM_P([1 : end 1], 2),...
        'r', 'LineWidth', 1)

    % Making the figure full screen programmatically
    set(gcf, 'Position', get(0, 'Screensize'));
                             
    % Releasing the figure that was held in function VelDomPlot
    hold off
end