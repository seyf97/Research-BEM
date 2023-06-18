function [mesh, assembly, values] = BEM_StokesPartDrag(mesh)
    % Traction and Drag Force on the Particle
    
    
    %%
    %{
    ***********************************************************************
        Assembly of System Matrices
    ***********************************************************************
    Generating the G and H matrices which will be repartitioned as a part 
    of the formulation process


        H * u = G * t


                 -  -
                | u0 |
            u = |    |
                | uP |
                 -  -

            u : (2 * N) x 1

            u0 : (2 * N0) x 1

            uP : (2 * NP) x 1


                 -  -
                | t0 |
            t = |    |
                | tP |
                 -  -

            t : (2 * N) x 1

            t0 : (2 * N0) x 1

            tP : (2 * NP) x 1   
    %}
    assembly = BEM_StokesPartAssem0(mesh);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nChannel Boundary Assembly done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till the 
    current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***********************************************************************
        Particle Boundary Assembly
    ***********************************************************************
    %}
    assembly = BEM_StokesPartAssemP(mesh, assembly);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nParticle Boundary Assembly done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till the 
    current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***************************************************************************
        Boudary Value Calculation
    ***************************************************************************

    We find boundary values of the channel only once and use them throughout
    the program because they do not chage significatly in case the particles 
    interact with the flow.


    %}
    
    % Initializing the values structure
    values = struct;
    
    %
    values = BEM_StokesPartBV(mesh, assembly, values);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nBoundary Value Calculation done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till 
    the current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***********************************************************************
     Plotting Velocity on Channel
    ***********************************************************************
    %}    
    BEM_StokesPartVel0Plot(mesh, values)
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nChannel Velocity Plot done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till 
    the current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***********************************************************************
       Re-Partition
    ***********************************************************************
    %}
    assembly = BEM_StokesPartPartition(mesh, assembly);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nRe-Partition done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till 
    the current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***********************************************************************
    Traction on the particle
    ***********************************************************************
    %}
    values.tP = assembly.BInv * assembly.C * values.u0;
    
    %{
    The forces from directions and the torque about the CG of 
    the particle on the particle is easily calculated from the 
    matrices at hand when the CG of the particle is at a given 
    pair of coordinates.


                 -  -
                | f1 |
           fB = | f2 |
                |  m |
                 -  -

    B : standing for body

    f1 : x component of the resultant force at the CG of the 
        particle

    f2 : y component of the resultant force at the CG of the 
        particle

    m : torque applied to the particle about its CG where the 
        positive direction is taken as counter-clockwise    
    
    We negate the force since it is positive when the stress 
    is directed inside the system. However, since we desire to see the 
    effect of the flow on the particle we are reversing the direction for 
    illustration.
    %}
    values.fB = - assembly.R * values.u0;
    
    %{
        Drag Wall Correction Factor
        (or Wall Correction Factor of the Drag Force
         or Nondimensional Drag Force)

    Adopted from [4]

             (Drag Force) / [(Maximum Fully Developed Laminar Speed) *
             (dynamic viscosity)]

    It is crucial to use the maximum value of the horizontal component
    of the fully developed parabolic velocity profile (rather than the 
    inlet velocity) in this formula in order to get consistent 
    non-dimensional force values with the reference.
    
    non : Nondimensional
    %}
    values.nonfBx = values.fB(1 : 3 : end) / (max(values.u0) * mesh.mu);
    
    
    %%
    %{
    ***********************************************************************
     Plotting Traction on Particle
    ***********************************************************************
    %}
    BEM_StokesPartTracP_Plot(mesh, values)
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nParticle Traction Plot done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till 
    the current section (inclusive)
    %}
    toc
    
    
    %%
    % Generating a domain mesh
    mesh = BEM_StokesPartMeshDom(mesh);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nDomain Mesh Generation done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till 
    the current section (inclusive)
    %}
    toc
    
    %%
    %{
    ***********************************************************************
        Velocity Field in the Computational Domain
    ***********************************************************************
    %}
    
    % Generating the velocity field
    [assembly, values] = BEM_StokesPartVelDom(mesh, assembly, values);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nVelocity Field Generation done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till 
    the current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***********************************************************************
        Plotting the Velocity Field in the Computational Domain
    ***********************************************************************
    %}
    BEM_StokesPartVelDomPlot(mesh, values);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nVelocity Field Plot done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till the 
    current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***********************************************************************
        Pressure Field in the Computational Domain
    ***********************************************************************
    %}
    
    % Generating the pressure field
    [assembly, values] = BEM_StokesPartPresDom(mesh, assembly, values);
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nPressure Field Generation done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till 
    the current section (inclusive)
    %}
    toc
    
    
    %%
    %{
    ***********************************************************************
        Plotting the Pressure Field in the Computational Domain
    ***********************************************************************
    %}
    mesh = BEM_StokesPartPresDomPlot(mesh, assembly, values);

    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nPressure Field Plot done...\n')

    %{
    Finalizing the command measuring the runtime of the program till the 
    current section (inclusive)
    %}
    toc
end