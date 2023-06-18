function mesh = BEM_StokesPartMeshPFFBounCheck(mesh)

%This function checks whether any point on the particle crosses the upper
%segment of the horizontal pinched part or the upper part of the inlet
%channel.

%Initially the line equation between points 9 and 10 must be determined,
%which is y1 = ( (x - x9)(y10-y9)/(x10-x9) ) + y9
%If the the y coordinate of a X1_P point lies below the y value of the
%above equation, the point is under the boundary and there is no problem.
%If the y coordinate is larger than the y of equation, the molecule is
%off-bounds and an error msg will be given.
    
%Counter is to check if all points has passed the pinched segment. If yes,
%there is no need to enter the loop
cnt = 0;

%If all the points has passed the pinched segment, the function will
%by-pass the for loop and no time will be wasted in this funct.
%Also if the channel is not PFF, the function will directly by-pass the
%loop
if mesh.isPffBoun && mesh.isPFF
    
%This for loop will be checking for every single point on the molecule
for i = 1 : mesh.N_E_P
    
    %Checking if a single point is behind point 9
    if mesh.X1_P(i,1) < mesh.X_Vc(9,1)
        
        %checking if the point is under the line y1 (must not be LARGER,
        %not even equal. If equal, on the boundary
        
        if ((mesh.X1_P(i,1) - mesh.X_Vc(9,1))*(mesh.X_Vc(10,2) - mesh.X_Vc(9,2))/(mesh.X_Vc(10,1) - mesh.X_Vc(9,1)) ) + ...
                mesh.X_Vc(9,2) <= mesh.X1_P(i,2)
            error('The particle is off-limits. Refine timeStep.'); 
        end
            
            %Checking if a point is between 8 and 9
    elseif mesh.X1_P(i,1) >= mesh.X_Vc(9,1) && mesh.X1_P(i,1) <= mesh.X_Vc(8,1)
            
        %checking if the point is under the the y value of y9 (or y8). If
        %it is equal or greater, error.
        if mesh.X_Vc(9,2) <= mesh.X1_P(i,2) 
            error('The particle is off-limits. Refine timeStep.'); 
        end
        
        %Checking if a point passed the x8 value. If all the points pass
        %this value, the function will bypass the for loop and this
        %function will no longer check for boundaries.
        
    elseif mesh.X1_P(i,1) > mesh.X_Vc(8,1)
        cnt = cnt + 1;
    end
end
end


%If all the points passed, the statement will be false and the funct will
%not enter the for loop
if cnt == mesh.N_E_P
    mesh.isPffBoun = false;
end

end

