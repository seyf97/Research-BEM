function mesh = BEM_StokesPartMeshDom(mesh)
    % Domain Mesh Generation
    
    
    %%
    %{
    ***********************************************************************
        Creating a Scalar PDE Model
    ***********************************************************************
    
    This model is going to be used to build the domain mesh into.
    %}
    mesh.model = createpde;
    
    
    %%
    %{
    ***********************************************************************
        Polygon Domain Mesh Generation Parameters
    ***********************************************************************
    
    We do not including the same point twice since the app connects the 
    first and last point automatically. In case, you include two same 
    points you will face the error message saying polygon points cannot 
    interact or overlap
    %}
    
    % x coordinates of the vertices of the polygon
    Polyx = mesh.X_Vc(:, 1);
    
    % y coordinates of the vertices of the polygon
    Polyy = mesh.X_Vc(:, 2);
    
    %{
        Matrix containing information about the polygon
    
    The first entry of the matrix, here 2, indicates to the function decsg
    which will be used to set up the geometry for mesh generation that the
    following coordinates belong to the vertices of a polygon.
    
    The second entry specifies the number of sides (line segments) the 
    polygon will have.
    
    The points go along the polygon in the counter-clockwise direction.
    %}
    Poly = [2; mesh.N_S_Poly; Polyx; Polyy];
    
    %{
    Name of the polygon for being placed in the name space and for use in
    set formula
    %}
    namePoly = 'Poly0';
    
    
    %%
    %{
    ***********************************************************************
        Circle Domain Mesh Generation Parameters
    ***********************************************************************
    %}
    
    %{
        Matrix containing information about the circle

    The first entry of the matrix, here 1, indicates to the function decsg
    which will be used to set up the geometry for mesh generation that the
    following coordinates and the radii belong to the circles.
    %}
    
    % Number of the circles present inside the domain
    mesh.N_Circ = mesh.N_P + mesh.N_P_Small + mesh.N_Post;
    
    %{
    Creating a matrix which stores the centers of all the cirles present
    inside the domain
    %}
    mesh.X_CG_O = [mesh.X_CG; mesh.X_CG_Small; mesh.X_CG_Post];
    
    %{
    Creating a matrix which stores the radii of all the cirles present
    inside the domain
    %}
    mesh.RadO = [mesh.Rad; mesh.RadSmall * ones(mesh.N_P_Small, 1);...
        mesh.postRad * ones(mesh.N_Post, 1)];
    
    % Initializing the matrix
    Circ = zeros(mesh.N_Circ, 4);
    
    % Adding information about each circle
    for i = 1 : mesh.N_Circ
        Circ(i, :) = [1 mesh.X_CG_O(i, 1) mesh.X_CG_O(i, 2) mesh.RadO(i)];
    end
    
    %{
    Transposing Circ for making it ready to be concatanated with Poly to
    form the geometric description
    %}
    Circ = Circ';
    
    %{
    Padding the matrix Circ with zeros since we need to combine it with
    the matrix Poly for forming the geoDescr matrix which will be input
    to the function decsg
    %}
    Circ = [Circ; zeros(length(Poly) - 4, mesh.N_Circ)];
    
    %{
    We initialize the name space for circles by prescribing the name of the
    first circle.
    %}
    mesh.nameCirc = ['Circ' num2str(1)];
    
    % We add the names of the circles to the name space
    for i = 1 : mesh.N_Circ - 1
        mesh.nameCirc = char(mesh.nameCirc, ['Circ' num2str(i + 1)]);
    end
    
    
    %%
    %{
    ***********************************************************************
    Forming Parameters for Geometry Decomposition
    ***********************************************************************
    
    We need 3 parameters in this part and they are created by using the
    previously formed parameters above. They are:
    
    1) Geometric Description of the Problem
    
    2) Name Space of the Geometry of the Problem
    
    3) Set Formula Defining the Geometry of the Problem
    %}
    
    % Geometric Description
    mesh.geoDescr = [Poly Circ];
    
    %{
    Name Space finalized by combining the names for the polygon and the
    circles
    %}
    mesh.nameSpace = char(namePoly, mesh.nameCirc);
    
    %{
    Adding the polygon and putting a minus for subtracting the
    upcoming circles which will be added together inside
    parantheses
    %}
    mesh.setFormula = [mesh.nameSpace(1, :) '-' '('];
    
    % Forming the set formula by utilizing the name space
    for i = 1 : mesh.N_Circ
        % Adding the circles
        mesh.setFormula = [mesh.setFormula mesh.nameSpace(i + 1, :) '+']; %#ok<*AGROW>
    end
    
    %  If all the circles are done, delete the trailing plus sign
    mesh.setFormula(end) = '';
    
    % Close the parantheses
    mesh.setFormula = [mesh.setFormula ')'];
    
    %{
    Making the names readable columnwise since that is how the function
    decsg processes
    %}
    mesh.nameSpace = mesh.nameSpace';
    
    
    %%
    %{
    ***********************************************************************
        Mesh Generation and Illustration
    ***********************************************************************
    %}
    
    %{
        Decomposing Constructive Solid Geometry
    
    dl : Decomposed Geometry Matrix containing information about the
    minimal regions
    
    bt : Boolean Table containing information about the borders between 
    minimal regions
    %}
    [dl, bt] = decsg(mesh.geoDescr, mesh.setFormula, mesh.nameSpace);
    
    %{
        Deleting borders between minimal regions
    
    For example, if there are line segments coinciding or intersecting the
    circles inside the polygon, csgdel deletes them.
    %}
    [dl, ~] = csgdel(dl, bt);
    
    % Including the geometry inside the model structure
    geometryFromEdges(mesh.model, dl);
    
    %{
    Generating the mesh for the geometry specified
    
    Information about the mesh such as the minimum and maximum element 
    sizes can be accessed from model.Mesh.
    %}
    generateMesh(mesh.model, 'Hmax', mesh.domMeshFrac * mesh.E_Size);
    
    %{
    Now that the mesh is stored inside the model structure, we access it 
    from there and store the nodes in a matrix to be output.
    
    Note that this matrix contains the nodes which coincide the boundary
    and since we need the nodes to be strictly inside the domain, we remove
    the nodes coinciding the boundary.
    %}
    mesh.X_D = mesh.model.Mesh.Nodes;
    
    %{
    Transposing the node matrix since the convention in the program for
    position matrices are columnwise 
    %}
    mesh.X_D = mesh.X_D';
    
    %{
    Finding the IDs of the nodes on the edges of the polygon and circles
    %}
    mesh.edgeNodeIDs = findNodes(mesh.model.Mesh, 'region', 'Edge',...
        1 : mesh.model.Geometry.NumEdges);
    
    % Forming the matrix containing the edge nodes
    mesh.X_D_Edge = mesh.model.Mesh.Nodes(:, mesh.edgeNodeIDs);
    
    %{
    Transposing the node matrix since the convention in the program for
    position matrices are columnwise 
    %}
    mesh.X_D_Edge = mesh.X_D_Edge';
    
    % Removing the nodes which coincide the boundary
    mesh.X_D = setdiff(mesh.X_D, mesh.X_D_Edge, 'rows');
    
    % Opening a new figure for plotting the domain mesh
    figure(mesh.figNumDomMesh)
    
    % Plotting the Mesh Structure
    pdeplot(mesh.model)
    
    % Holding the current figure
    hold on
    
    % Plotting the nodes inside the domain
    scatter(mesh.X_D(:, 1), mesh.X_D(:, 2))

    % Making the figure full screen programmatically
    set(gcf, 'Position', get(0, 'Screensize'));
    
    % Releasing the current figure
    hold off
    
    % Number of points/nodes generated in the domain
    mesh.N_D = length(mesh.X_D);
end