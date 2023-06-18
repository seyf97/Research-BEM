function BEM_StokesPartTrajMov(mesh, values)
    % Particle Trajectory Movie
    %
    %

    % Opening the appropriate figure
    figure(mesh.figNumTrajMov)

    % Plotting the channel
    BEM_StokesPartMeshBounPlot0(mesh);

    % Deleting the potential trailing zeros
    values.allX_CG = values.allX_CG(:, :, 1 : find(values.allX_CG(1, 1, :), 1, 'last'));

    % Creating the movie and rewinding it the prescribed times
    for i = 1 : size(values.allX_CG, 1)
        comet(values.allX_CG(i, 1, :), values.allX_CG(i, 2, :))
    end
end