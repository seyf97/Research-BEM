function mesh = BEM_StokesPartMesh0(mesh)    
    % Channel Boundary Mesh Generation
    % Forming the parameters related to the channel that are necessary for 
    % the assembly of the system matrices
    
    
    % Discretization in case of Simple Rectangular Channel
    if mesh.isSRC
        mesh = BEM_StokesPartMeshSRC(mesh);
    
    % Discretization in case of Hydrodynamic Particle Separation
    elseif mesh.isHPS %#ok<*UNRCH>
        mesh = BEM_StokesPartMeshHPS(mesh);
    
    % Discretization in case of Pinched Flow Fractionation
    elseif mesh.isPFF
        mesh = BEM_StokesPartMeshPFF(mesh);
    end
    
    % Discretization in case of Deterministic Lateral Displacement
    if mesh.isDLD
        mesh = BEM_StokesPartMeshDLD(mesh);
    end
end