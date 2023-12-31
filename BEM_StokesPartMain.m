%{
***************************************************************************
    Preface
***************************************************************************

---------------------------------------------------------------------------
  Coding Style

Each comment addresses the line of code below it.

Each section of code is separated from another by 2 new lines.

Inside every section a line of code is separated by a single new line
from another line of code.

Each section has a title contained within 2 lines of star symbols.

Subsections are contained within 2 lines of hyphen symbol.

If there is a single heading/title for the code of line coming after the
comment inside a section, do not indent the heading. If other
explanations are included in the comment, indent the heading for the 
purposes of visual emphasis.

Formulas are indented three times inside a comment.

If a variable is adopted from a source, indicate the reference and add it
to the end of the program in the References section for later access.
---------------------------------------------------------------------------


---------------------------------------------------------------------------
  Mathematical Notions

A vector refers to an array of size N x 1 where N is a natural number.

A matrix is an array of size M x N where M & N are both natural numbers.
---------------------------------------------------------------------------


---------------------------------------------------------------------------
  Physical Notions

Channel : static set-up that the fluid flows through which also includes
all the contractions and obstacles staticly present inside the set-up
(designated by 0 throughout the program)

Particle : Moving or stationary entities inside the channel (designated by
p or P throughout the program)
---------------------------------------------------------------------------
%}


%%
%{
***************************************************************************
    Preparing Workstations
***************************************************************************
%}

% Clearing the Variable Workspace
clear

% Clearing the Command Window
clc

% Closing all the possibly open figures
close all

% Initializing the command to measure the runtime of the program
tic


%%
%{
***************************************************************************
    Find Drag Force
***************************************************************************

    true : find the force acting on the particle

    false : move the particle with the flow

Note that this program finds the force by exploiting the fact that the
particle is stationary at the initial center of gravity given and that the
program terminates right after finding the force.
%}
mesh.findForce = ~true;


%%
%{
***************************************************************************
    Figures
***************************************************************************
%}

%{
Figure on which the state of the boundaries of the microchannel will be 
illustrated
%}
mesh.figNumBounMesh = 1;

%{
Figure on which the state of the domain of the microchannel will be 
illustrated
%}
mesh.figNumDomMesh = 2;

% Figure demonstrating the velocity field in the domain
mesh.figNumVelField = 3;

% Figure demonstrating the pressure field in the domain
mesh.figNumPresField = 4;

% Figure in which the trajectory movie will be illustrated
mesh.figNumTrajMov = 5;


% In case particulate flow will be observed
%%
if ~mesh.findForce
    %{
    ***********************************************************************
        Simulation Time Parameters
    ***********************************************************************
    %}

    %{
        Time Step [s]

    Time step is going to be used to update the position of the particle 
    once its rigid-body motion parameters (the translational velocity 
    components, x and y, and the rotational/angular velocity) are known.

    Throughout a time step the rigid-body motion parameters of the particle
    are assumed to be constant.
    %}
    mesh.timeStep = .01;
    
    % Duration of the Simulation [s]
    mesh.simTime = 1.5;
    
    %Counter to generate random molecules only once (Do not adjust)
    mesh.counter = 0;

    %{
    Time variable that will be updated and compared with simTime at every
    movement of the particle(s)
    %}
    mesh.time = mesh.timeStep;
end
%%
%{
***************************************************************************
SEYFO UPDATES AS OF 03/12/18

UPDATES:

1)RANDOM MOL GENERATOR ADDED
1a)THE XCG VALUES ARE ALL BASED ON THE REAL XCG OF THE RANDOM MOLECULE. 

2)ELEMENT NUMBERS WILL BE ASSIGNED AUTOMATICALLY BY A FRACTION OF E_LENGTH

3)IF MOLECULE TOUCHES OR PASSES BOUNDARY WHILE IN THE PINCHED SEGMENT,
ERROR MESSAGE WILL BE GIVEN.
%}

%Statement to check whether the molecule has passed the pinched segment. If
%it is true, it is inside the pinched segment of the channel. If it passes
%the pinched segment, it will become false.
mesh.isPffBoun = true;



%%
%{
***************************************************************************
    Fluid Properties
***************************************************************************
%}

% Dynamic viscosity of the fluid [Pa * s]
mesh.mu = 10 ^ - 3;

% Density of the fluid [kg / m ^ 3]
mesh.rho = 10 ^ 3;

% Kinematic viscosity of the fluid [m ^ 2 / s]
mesh.nu = mesh.mu / mesh.rho;
    
% Horizontal Component of Inlet Velocity [m / s]
mesh.vIn = 40 * 10 ^ - 6;

% Horizontal Component of Outlet Traction [N / m ^ 2]
mesh.tOut = 0;


%%
%{
***************************************************************************
    Channel Boundary Discretization and Geometry Specifications
***************************************************************************

Only one of the boolean variables in this section can be true at a time
with the exception of the pair containing isDLD.
%}

%{
    Is there a Simple Rectangular Channel?

SRC : Simple Rectangular Channel
%}
mesh.isSRC = true;

%{
    Is there a structure for Hydrodynamic Particle Separation?

HPS : Hydrodynamic Particle Separation
%}
mesh.isHPS = ~true;

%{
    Is there a structure for Pinched Flow Fractionation?

PFF : Pinched Flow Fractionation
%}
mesh.isPFF = ~true;

%{
    Is there a structure for Deterministic Lateral Displacement?

DLD : Deterministic Lateral Displacement
%}
mesh.isDLD = ~true;

% Element Size
mesh.E_Size = 10 * 10 ^ - 6;

%The size of the circular elemetns will be E_L/ECMultiplier
mesh.E_C_Multiplier = 20;

% What factor of boundary element size should the domain mesh elements be?
mesh.domMeshFrac = 2;

% Parameters for Simple Rectangular Channel
if mesh.isSRC
    %{
        Length of the bottom or top sides of the channel wall [m]

    In case of a contraction placed for Hydrodynamic Separation simulation,
    this length is partitioned inside the Mesh code to 5 pieces.

    L : Length
    %}
    mesh.Lx = 180 * 10 ^ - 6;
    
    % Length of the left or right sides of the channel wall [m]
    mesh.Ly = 50 * 10 ^ - 6;
    
% Parameters for Hydrodynamic Particle Separation
elseif mesh.isHPS %#ok<*UNRCH>
    % Length of the bottom or top sides of the channel wall [m]
    mesh.Lx = 500 * 10 ^ - 6;
    
    % Length of the left or right sides of the channel wall [m]
    mesh.Ly = 100 * 10 ^ - 6;
    
    % Contraction Roof Length [m]
    mesh.L_Gapx = 30 * 10 ^ - 6;
    
    % Contraction Gap Length [m]
    mesh.L_Gapy = 30 * 10 ^ - 6;
    
    % Contraction Base Length [m]
    mesh.L_Base = 35 * 10 ^ - 6;
    
% Parameters for Pinched Flow Fractionation
elseif mesh.isPFF
    % Horizontal length of the broadened segment of the channel
    mesh.L_Brox = 2000 * 10 ^ - 6;
    
    % Vertical length of the broadened segment of the channel
    mesh.L_Broy = 1000 * 10 ^ - 6;
    
    % Horizontal length of the pinched segment of the channel
    mesh.L_Pinx = 100 * 10 ^ - 6;
    
    % Vertical length of the pinched segment of the channel
    mesh.L_Piny = 50 * 10 ^ - 6;
    
    % Vertical length of both inlets to the channel
    mesh.L_Iny = 125 * 10 ^ - 6;
    
    %{
    Angle of the diffuser part (transition from pinched segment to 
    broadened segment) of the PFF channel
    %}
    mesh.diffuserAng = pi;
    
    %{
    Angle of the nozzle part (transition from inlet to pinched segment) of 
    the PFF channel
    %}
    mesh.nozzleAng = pi / 3;
    
    % Depth of the channel
    mesh.Lz = 50 * 10 ^ - 6;
    
    % Volumetric Flow Rate of the pure fluid [L / h]
    mesh.Q_In0 = 120 * 10 ^ - 6;
    
    % Volumetric Flow Rate of the fluid containing particles [L / h]
    mesh.Q_InP = 20 * 10 ^ - 6;
    
    % Volumetric Flow Rate of the pure fluid [m ^ 3 / s]
    mesh.Q_In0 = mesh.Q_In0 * 10 ^ - 3 / 3600;
    
    % Volumetric Flow Rate of the fluid containing particles [m ^ 3 / s]
    mesh.Q_InP = mesh.Q_InP * 10 ^ - 3 / 3600;
    
    % Inlet Velocity of the pure fluid [m / s]
    mesh.vIn0 = mesh.Q_In0 / (mesh.L_Iny * mesh.Lz);

    % Inlet Velocity of the fluid containing particles [m / s]
    mesh.vInP = mesh.Q_InP / (mesh.L_Iny * mesh.Lz);

    % Horizontal Component of Outlet Traction
    mesh.tOut = 0;
end

% Parameters for Deterministic Lateral Displacement
if mesh.isDLD
    % Column-to-column Spacing
    mesh.offsetx = 8 * 10 ^ - 6;
    
    % Vertical Obstacle Spacing
    mesh.offsety = 8 * 10 ^ - 6;
    
    % Gap Width
    mesh.gapy = 1.6 * 10 ^ - 6;
    
    %{
        Post Radius

    Rad : Radius
    %}
    mesh.postRad = (mesh.offsety - mesh.gapy) / 2;
    
    %{
        Shift Fraction

    Frac : Fraction
    %}
    mesh.shiftFrac = .1;
    
    % Vertical distance each column is shifted from the previous column
    mesh.shifty = mesh.shiftFrac * mesh.offsety;
    
    %{
        Post Array Length in x direction

    Arr : Array
    %}
    mesh.postArrLx = mesh.Lx;
    
    % Post Array Length in y direction
    mesh.postArrLy = mesh.Ly;
    
    %{
    Coordinates of the 1st post at the bottom left corner of the post array
    %}
    mesh.X_CG_1 = [mesh.postArrLx * 20 / 180 mesh.postArrLy * 5 / 50];
else
    
    %{
        Post Radius

    Rad : Radius
    %}
    mesh.postRad = 0;
    
    %
    mesh.N_Post = 0;
    
    %
    mesh.X_CG_Post = [];
end 


%%
%{
***************************************************************************
    Particle Boundary Discretization and Geometry Specifications
***************************************************************************
%}

%{
    Number of Particles

P : Particle

N_P : 1 x 1
%}
mesh.N_P = 1;

% In case of PFF
if mesh.isPFF
    % Radii of the particles [m]
    mesh.Rad = 10 * 10 ^ - 6;
    
    % Initial Positions of the centers of gravity of the particles
    mesh.X_CG = [4 / 5 * (mesh.L_Broy / 2 - mesh.L_Iny) /...
        tan(mesh.nozzleAng / 2)...
        4 / 5 * (mesh.L_Iny - mesh.L_Broy / 2) + ...
         (mesh.L_Broy - mesh.L_Iny / 2)];
    
    % Initial angles the particles make with the x axis [Radians]
    mesh.rotAng = zeros(mesh.N_P, 1);
    
% In case of HPS
elseif mesh.isHPS
    
    %{
    Radius of the particle

    Rad : Radius

    Rad : N_P x 1
    %}
    mesh.Rad = 5 * 10 ^ - 6;
    
    %{
    Initial Position of the center of gravity of the particle

    CG : Center of Gravity

    X_CG : N_P x 2 matrix
    %}
    mesh.X_CG = [mesh.Lx * 15 / 50 mesh.Ly * 1.5 / 5];
    
    %{
    Initial angle the particle makes with the x axis

    rot : Rotation

    Ang : Angle

    rotAng : N_P x 1
    %}
    mesh.rotAng = 0;
    
% In case of the flow of a single particle
elseif mesh.isSRC
    
    %{
    Radius of the particle

    Rad : Radius

    Rad : N_P x 1
    %}
    mesh.Rad = 6.95 * 10 ^ - 6;
    
    %{
    Initial Position of the center of gravity of the particle

    CG : Center of Gravity

    X_CG : N_P x 2 matrix
    %}
    mesh.X_CG = [mesh.Lx /2 mesh.Ly/4];
    
    %{
    Initial angle the particle makes with the x axis

    rot : Rotation

    Ang : Angle

    rotAng : N_P x 1
    %}
    mesh.rotAng = 0;
    
% In case of 8 particles
else
    % Radii of the particles [m]
    mesh.Rad = 10 * [.40 / 2 * 10 ^ - 6;...
                .50 / 2 * 10 ^ - 6;...
                .60 / 2 * 10 ^ - 6;...
                .70 / 2 * 10 ^ - 6;...
                .80 / 2 * 10 ^ - 6;...
                .90 / 2 * 10 ^ - 6;...
               1.03 / 2 * 10 ^ - 6;...
               1.10 / 2 * 10 ^ - 6];
    
    % Initial Positions of the centers of gravity of the particles
    mesh.X_CG = [mesh.Lx * 190 / 500 * ones(mesh.N_P, 1)...
                 mesh.Ly * (10 : 25 : 190)' / 200];
    
    % Initial angles the particles make with the x axis [Radians]
    mesh.rotAng = zeros(mesh.N_P, 1);
end

%{
    Number of elements on the particles

We multiply whatever comes out of the ceil function by 3 because
constructing a circle with less than 3 elements is not realistic and unique
%}
% mesh.N_E_P = 20 * ceil(2 * pi * mesh.Rad / mesh.E_Size);

%{
    Attached Particles

Mol : Molecule
%}

% Are there attached uniformly distributed particles?
mesh.isMolUni = false;

%Are there attached randomly distributed particles?
mesh.isMolRand = false;

% In case the particles are molecules
if mesh.isMolUni || mesh.isMolRand
    % Number of attached particles
    mesh.N_P_Small = 12;
    
    %{
    Number of elements on the small attached particles
    
    Note that all the attached particles will have the same number of
    elements on them and that this number indicates the number of elements 
    on a single small particle. Thus, overall there will be (N_P_Small *
    N_E_P_Small) elements on the small particles.
    %}
% % % % % % % % % % % %     mesh.N_E_P_Small = 16 * mesh.N_P_Small;
    
    % Number of elements on the large particle
% % % % % % % % % %     mesh.N_E_P_Large = 64 * mesh.N_P_Small;
    
    %
% % % % % % % % % % % % %     mesh.N_E_P = mesh.N_E_P_Large + mesh.N_E_P_Small * mesh.N_P_Small;
    
    %{
    Radius of the attached particles
    
    Note that all the attached particles will have the same radius.
    %}
    mesh.RadSmall = 1 * 10 ^ - 6;
    
    % Penetration amount [m]
    mesh.penAm = .05 * 10 ^ - 6;
else
    
    % Number of attached particles
    mesh.N_P_Small = 0;
    
    %{
    Radius of the attached particles
    
    Note that all the attached particles will have the same radius.
    %}
    mesh.RadSmall = 0;
    
    %
    mesh.X_CG_Small = [];
end

%CREATING THE NUMBER OF ELEMENTS FOR A SINGLE MOLECULE ACCORDING TO THE
%INITAL ELEMENT LENGTH

mesh = BEM_StokesPartMeshPElNum(mesh);
%%
%{
***************************************************************************
    Fluid Properties and Flow Conditions
***************************************************************************
%}

if mesh.isSRC
    %{
        Characteristic dimensions of the particles

    Adopted from [4]

    a : N_P x 1
    %}
    mesh.a = mesh.Rad;
    
    %{
        Characteristic width (half-width) of the channel

    Adopted from [4]
    %}
    mesh.b = mesh.Ly / 2;
    
    %{
        Characteristic Aspect Ratios for the channel and the particles

    Adopted from [4]

            (Radius of the particle) / (Half-width of the channel)

    k : N_P x 1
    %}
    values.k = mesh.a / mesh.b;
    
    %{
        Reynolds number characterizing the flow

    Adopted from [4]

            (Horizontal Component of Inlet Velocity) * 
            * (Characteristic Particle Dimension) / (Kinematic Viscosity)

    Re : Reynolds

    Re : N_P x 1
    %}
    values.Re = mesh.vIn * mesh.a / mesh.nu;
end


%%
%{
***************************************************************************
   Channel Boundary Mesh Generation
***************************************************************************
%}
mesh = BEM_StokesPartMesh0(mesh);

% Prompting the user about the completion of the operation in the section
fprintf('\nChannel Boundary Mesh Generation done...\n')

%{
Finalizing the command measuring the runtime of the program till the 
current section (inclusive)
%}
toc


%%
%{
***************************************************************************
   Particle Boundary Mesh Generation
***************************************************************************

Forming the parameters related to the particle that are necessary for the 
assembly of the system matrices
%}
%Making sure this first mesh is the first iteration
%values a is just to make sure values structure is created. just ignore it lol
values.a = 0;
mesh.iter = false;
mesh = BEM_StokesPartMeshP(mesh,values);
mesh.iter = true;


% Total number of elements
mesh.N_E = mesh.N_E_0 + mesh.N_E_P_O;

% Prompting the user about the completion of the operation in the section
fprintf('\nParticle Boundary Mesh Generation done...\n')

%{
Finalizing the command measuring the runtime of the program till the 
current section (inclusive)
%}
toc


%%
%{
***************************************************************************
   Particle Boundary Mesh Plot
***************************************************************************

Plotting the elements on particle boundary
%}
BEM_StokesPartMeshBounPlot(mesh);


%Checking if the particle is inside the boundary limits
mesh = BEM_StokesPartMeshPFFBounCheck(mesh);

% Prompting the user about the completion of the operation in the section
fprintf('\nBoundary Mesh Plot done...\n')

%{
Finalizing the command measuring the runtime of the program till the 
current section (inclusive)
%}
toc


%%
% In case particulate flow will be observed
if ~mesh.findForce
    %%
    %{
    ***********************************************************************
        Update Storage Matrices
    ***********************************************************************
    %}
    
    %{
        Initializing the matrix containing the positions of the center of 
    gravity of the particle(s)

    allX_CG : Np x 2 x (simTime / timeStep)
    %}
    values.allX_CG = zeros(mesh.N_P, 2, mesh.simTime / mesh.timeStep);
    
    % Setting the first positions to the prescribed positions
    values.allX_CG(:, :, 1) = mesh.X_CG;
    
    %{
        Initializing the matrix containing the rotation angles of the 
    particles

    allRotAngs : Np x (simTime / timeStep)
    %}
    values.allRotAngs = zeros(mesh.N_P, mesh.simTime / mesh.timeStep);
    
    % Setting the first rotation angles to the prescribed angles
    values.allRotAngs(:, 1) = mesh.rotAng;
    
    %{
        Initializing the matrix for storing rigid-body motion parameters 
    of the particles

    allUB : 3 * Np x (simTime / timeStep)
    %}
    values.allUB = zeros(3 * mesh.N_P, mesh.simTime / mesh.timeStep);
end

%%
% In case Drag Force on the stationary particle is sought
if mesh.findForce
    [mesh, assembly, values] = BEM_StokesPartDrag(mesh);
    
% In case particulate flow will be observed
else
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
    As long as the simulation time is not surpassed and the trajectories of the
    particles are confined within the channel boundaries, the program keeps 
    updating the position of the particle.
    %}
    while mesh.time < mesh.simTime && all(all(max(mesh.X_Vc) > mesh.X1_P))
        %%
        %{
        *******************************************************************
            Particle Boundary Assembly
        *******************************************************************
        %}
        assembly = BEM_StokesPartAssemP(mesh, assembly);
        
        %{
        Prompting the user about the completion of the operation in the 
        section
        %}
        fprintf('\nParticle Boundary Assembly done...\n')
        
        %{
        Finalizing the command measuring the runtime of the program till 
        the current section (inclusive)
        %}
        toc
        
        
        %%
        %{
        *******************************************************************
            Boudary Value Calculation
        *******************************************************************

        We find boundary values of the channel only once and use them 
        throughout the program because they do not chage significatly in 
        case the particles interact with the flow.


        %}
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
        *******************************************************************
         Plotting Velocity on Channel
        *******************************************************************
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
        Prompting the user about the completion of the operation in the section
        %}
        fprintf('\nRe-Partition done...\n')
        
        %{ 
        Finalizing the command measuring the runtime of the program till the
        current section (inclusive) 
        %}
        toc
        
        
        %%
        %{
        *******************************************************************
        Calculation of Rigid-body Motion Parameters
        *******************************************************************

        Solving the linear system to get the unknown values
        %}
        values = BEM_StokesPartRigBodyPars(assembly, values);
        
        % Adding the motion parameters to the storing matrix
        values.allUB(:, round(mesh.time / mesh.timeStep)) = values.uB;
        
        %{
        Prompting the user about the completion of the operation in the 
        section
        %}
        fprintf('\nCalculation of Rigid-body Motion Parameters done...\n')

        %{ 
        Finalizing the command measuring the runtime of the program till 
        the current section (inclusive) 
        %}
        toc
        
        
% % %         %%
% % %         %{
% % %         *******************************************************************
% % %             Domain Mesh Generation
% % %         *******************************************************************
% % %         %}
% % % 
% % %         % Generating a domain mesh
% % %         mesh = BEM_StokesPartMeshDom(mesh);
% % % 
% % %         %{
% % %         Prompting the user about the completion of the operation in the 
% % %         section
% % %         %}
% % %         fprintf('\nDomain Mesh Generation done...\n')
% % % 
% % %         %{
% % %         Finalizing the command measuring the runtime of the program till 
% % %         the current section (inclusive)
% % %         %}
% % %         toc
% % %         
% % %         %%
% % %         %{
% % %         *******************************************************************
% % %             Velocity Field in the Computational Domain
% % %         *******************************************************************
% % %         %}
% % % 
% % %         % Generating the velocity field
% % %         [assembly, values] = BEM_StokesPartVelDom(mesh, assembly, values);
% % %         
% % %         %{
% % %         Prompting the user about the completion of the operation in the 
% % %         section
% % %         %}
% % %         fprintf('\nVelocity Field Generation done...\n')
% % % 
% % %         %{
% % %         Finalizing the command measuring the runtime of the program till 
% % %         the current section (inclusive)
% % %         %}
% % %         toc
% % % 
% % % 
% % %         %%
% % %         %{
% % %         *******************************************************************
% % %             Plotting the Velocity Field in the Computational Domain
% % %         *******************************************************************
% % %         %}
% % %         BEM_StokesPartVelDomPlot(mesh, values);
% % % 
% % %         %{
% % %         Prompting the user about the completion of the operation in the 
% % %         section
% % %         %}
% % %         fprintf('\nVelocity Field Plot done...\n')
% % % 
% % %         %{
% % %         Finalizing the command measuring the runtime of the program till 
% % %         the current section (inclusive)
% % %         %}
% % %         toc
% % % 
% % % 
% % %         %%
% % %         %{
% % %         *******************************************************************
% % %             Pressure Field in the Computational Domain
% % %         *******************************************************************
% % %         %}
% % % 
% % %         % Generating the pressure field
% % %         [assembly, values] = BEM_StokesPartPresDom(mesh, assembly, values);
% % %         
% % %         %{
% % %         Prompting the user about the completion of the operation in the 
% % %         section
% % %         %}
% % %         fprintf('\nPressure Field Generation done...\n')
% % % 
% % %         %{
% % %         Finalizing the command measuring the runtime of the program till 
% % %         the current section (inclusive)
% % %         %}
% % %         toc
% % %         
% % %         
% % %         %%
% % %         %{
% % %         *******************************************************************
% % %             Plotting the Pressure Field in the Computational Domain
% % %         *******************************************************************
% % %         %}
% % %         mesh = BEM_StokesPartPresDomPlot(mesh, assembly, values);
% % %         
% % %         %{
% % %         Prompting the user about the completion of the operation in the 
% % %         section
% % %         %}
% % %         fprintf('\nPressure Field Plot done...\n')
% % %         
% % %         %{
% % %         Finalizing the command measuring the runtime of the program till 
% % %         the current section (inclusive)
% % %         %}
% % %         toc
% % %         
% % %         
        %%
        %{
        *******************************************************************
            Printing the Update Number
        *******************************************************************

        This number designates how many times the particle appeared in
        different spots inside the channel.
        %}
        fprintf('\nUpdate %d / %d \n', round(mesh.time / mesh.timeStep), round(mesh.simTime/mesh.timeStep))
        
% % %         %
% % %         drawnow
% % %         
% % %         % 
% % %         values.FrameArr1(round(mesh.time / mesh.timeStep)) =...
% % %             getframe(figure(1));
% % %         
% % %         %
% % %         values.FrameArr2(round(mesh.time / mesh.timeStep)) =...
% % %             getframe(figure(2));
% % %         
% % %         %
% % %         values.FrameArr3(round(mesh.time / mesh.timeStep)) =...
% % %             getframe(figure(3));
% % %         
% % %         %
% % %         values.FrameArr4(round(mesh.time / mesh.timeStep)) =...
% % %             getframe(figure(4));
        
        
        %%
        %{
        *******************************************************************
           Particle Position and Rotation Angle Update
        *******************************************************************
        %}
        
        %{
        Updating the x component of the position of the center of gravity 
        of the particles since we know the translational velocity now
        %}
        mesh.X_CG(:, 1) = mesh.X_CG(:, 1) + values.uB(1 : 3 : end) *...
            mesh.timeStep;
        
        %{
        Updating the y component of the position of the center of gravity 
        of the particles
        %}
        mesh.X_CG(:, 2) = mesh.X_CG(:, 2) + values.uB(2 : 3 : end) *...
            mesh.timeStep;
        
        
        % Storing the new position
        values.allX_CG(:, :, round(mesh.time / mesh.timeStep) + 1) =...
            mesh.X_CG;

        %{
        Updating the rotation angle since we know the rotational velocity
        now
        %}
        mesh.rotAng = mesh.rotAng + values.uB(3 : 3 : end) * mesh.timeStep;
        
        % Storing the new rotation angle
        values.allRotAngs(:, round(mesh.time / mesh.timeStep) + 1) =...
            mesh.rotAng;

        %{
        Prompting the user about the completion of the operation in the 
        section
        %}
        fprintf('\nParticle Position and Rotation Angle Update done...\n')

        %{ 
        Finalizing the command measuring the runtime of the program till 
        the current section (inclusive) 
        %}
        toc
        
        
        %%
        %{
        *******************************************************************
        Time Incrementation
        *******************************************************************
        %}
        mesh.time = mesh.time + mesh.timeStep;
        
        %
        values.runTime = mesh.simTime / mesh.time * toc;
        
        %
        values.runTimeLeft = values.runTime - toc;
        
        %
        fprintf('\n\t%4.2f%% done... \n \t%.0f Minutes to Go...\n',...
            toc / values.runTime * 100, values.runTimeLeft / 60)
        
        
        %%
        %{
        *******************************************************************
           Particle Boundary Mesh Generation
        *******************************************************************

        Forming the parameters related to the particle that are necessary 
        for the assembly of the system matrices
        %}
        mesh = BEM_StokesPartMeshP(mesh,values);

        %{
        Prompting the user about the completion of the operation in the 
        section
        %}
        fprintf('\nParticle Boundary Mesh Generation done...\n')

        %{
        Finalizing the command measuring the runtime of the program till 
        the current section (inclusive)
        %}
        toc

        
        %%
        %{
        *******************************************************************
           Particle Boundary Mesh Plot
        *******************************************************************

        Plotting the elements on particle boundary
        %}
        BEM_StokesPartMeshBounPlot(mesh);

        %{
        Prompting the user about the completion of the operation in the 
        section
        %}
        
        %Checking if the particle is inside limits
        mesh = BEM_StokesPartMeshPFFBounCheck(mesh);
        
        fprintf('\nBoundary Mesh Plot done...\n')

        %{
        Finalizing the command measuring the runtime of the program till 
        the current section (inclusive)
        %}
        toc
    end
    
    
    %%
    %{
    ***********************************************************************
        Movies
    ***********************************************************************
    %}
    
    %{
        Plotting the particle trajectory as a movie

    Traj : Trajectory

    Mov : Movie
    %}
    BEM_StokesPartTrajMov(mesh, values);
    
% % %     % Boun Mesh Movie
% % %     videoWriterFile = VideoWriter('Boundary Mesh Movie.avi');
% % %     
% % %     open(videoWriterFile)
% % %     
% % %     writeVideo(videoWriterFile, values.FrameArr1)
% % %     
% % %     close(videoWriterFile)
% % %     
% % %     % Domain Mesh Movie
% % %     videoWriterFile = VideoWriter('Domain Mesh Movie.avi');
% % %     
% % %     open(videoWriterFile)
% % %     
% % %     writeVideo(videoWriterFile, values.FrameArr2)
% % %     
% % %     close(videoWriterFile)
% % %     
% % %     % Velocity Field Movie
% % %     videoWriterFile = VideoWriter('Velocity Field Movie.avi');
% % %     
% % %     open(videoWriterFile)
% % %     
% % %     writeVideo(videoWriterFile, values.FrameArr3)
% % %     
% % %     close(videoWriterFile)
% % %     
% % %     % Pressure Field Movie
% % %     videoWriterFile = VideoWriter('Pressure Field Movie.avi');
% % %     
% % %     open(videoWriterFile)
% % %     
% % %     writeVideo(videoWriterFile, values.FrameArr4)
% % %     
% % %     close(videoWriterFile)
    
    %{
    Prompting the user about the completion of the operation in the 
    section
    %}
    fprintf('\nMovies done...\n')
    
    %{
    Finalizing the command measuring the runtime of the program till the 
    current section (inclusive)
    %}
    toc
end


%%
%{
***************************************************************************
    References
***************************************************************************

[1]: Kogl

[2]: Besim's Notes

[3]: Wrobel

[4]: Baranoglu, B.; Cetin, B. (2014): A particle flow specific boundary 
    element formulation for microfluidic applications. In Karayiannis, T.;
    Konig, C. S.; Balabani, S.(Eds): 4th Micro and Nano Flow Conference, 
    6�10 September 2014, no. 206. Brunel University.

[5]: A. Ben Richou, A. Ambaria, M. Lebeyc, and J.K. Nacirid. Drag force on 
    a circular cylinder mid- way between two parallel plates at very low 
    Reynolds numbers�Part 1: Poiseuille flow (numerical). Chemical 
    Engineering Science, 59(15): 3215�3222, 2004.

[6]: https://www.dextran.com/about-dextran/dextran-chemistry/
     physical-properties
%}


%%
% Prompting the user about the completion of the simulation
fprintf('\nSimulation over...\n')

% Finalizing the command measuring the runtime of the whole program
toc