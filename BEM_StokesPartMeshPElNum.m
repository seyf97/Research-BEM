function mesh = BEM_StokesPartMeshPElNum(mesh)

elSize = mesh.E_Size/mesh.E_C_Multiplier;

if mesh.isMolUni
    
    %!!! Only unknowns are num of elements on a single gap and small
    %circle!!!
    
    %Note that this angle will be the same for all the attached particles.
    mesh.angSmall = 2 * acos((mesh.Rad ^ 2 + (mesh.Rad + mesh.RadSmall -...
        mesh.penAm) ^ 2 - mesh.RadSmall ^ 2) / (2 * mesh.Rad * ...
        (mesh.Rad + mesh.RadSmall - mesh.penAm)));
    
    %Lambda=Angle btw the centers of two small circles
    angCG_Small=2*pi/mesh.N_P_Small;
    
    %Half the distance between the two coincident points of the Small cirlce
    %(S) and the Large circle (L)
    b=mesh.Rad*real(sin(mesh.angSmall));

    %Distance between the center of a small circle and Origin
    L=mesh.Rad*real(cos(mesh.angSmall)) + sqrt(mesh.RadSmall.^2 - b.^2);
    
    %Dist btw outer points of the L and S circle
    SOuterPt=L-mesh.RadSmall;
    LOuterPt=mesh.Rad;
    mesh.Dist=LOuterPt-SOuterPt;
    
    %Beta, the angle of the starting point of meshing of the small circle
    Beta=atan(sqrt(mesh.RadSmall.^2 - b.^2)/b);
    
    
    %**************** NUMBER OF ELEMENTS ON A SINGLE SMALL CIRCLE*****************
    mesh.N_E_P_Small = ceil((pi + 2*Beta)/(2*asin(elSize/(2*mesh.RadSmall))));
    
    
    %The angle of a single gap which will be meshed
    MeshAngle=angCG_Small-2*mesh.angSmall;
    
    %*****************NUMBER OF ELEMENTS ON A SINGLE GAP**********************
    GapNumEl = ceil(MeshAngle/(2*asin(elSize/(2*mesh.Rad))));
    
    
    %******************* NUMBER OF ELEMENTS ON ALL THE GAPS ***********************
    mesh.N_E_P_Large = GapNumEl*mesh.N_P_Small;
    
    
    %IN THE PART BELOW, IT IS ASSUMED THAT THERE IS ONLY A SINGLE MOLECULE
    %IN THE WHOLE SIMULATION
    %NUMBER OF ELEMENTS
    mesh.N_E_P_O = mesh.N_E_P_Large + mesh.N_P_Small*mesh.N_E_P_Small;
    mesh.N_E_P = mesh.N_E_P_O;

    %Matrix containing the indices of the last elements of each circle
    mesh.lastE_P = mesh.N_E_P;

    % Matrix containing the indices of the first elements of each circle
    mesh.firstE_P = 1;
    
    
    %IF IT IS ONLY A SMOOTH CIRCLE
elseif ~mesh.isMolRand && ~mesh.isMolUni
    
    mesh.N_E_P = ceil(pi/(asin(elSize/(2*mesh.Rad))));
    
    
    %IN THE PART BELOW, IT IS ASSUMED THAT THERE IS ONLY A SINGLE MOLECULE
    %IN THE WHOLE SIMULATION
    %NUMBER OF ELEMENTS
    mesh.N_E_P_O = mesh.N_E_P;

    %Matrix containing the indices of the last elements of each circle
    mesh.lastE_P = mesh.N_E_P;

    % Matrix containing the indices of the first elements of each circle
    mesh.firstE_P = 1;
    
end

end

    
  