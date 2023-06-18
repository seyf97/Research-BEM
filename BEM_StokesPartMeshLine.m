function [X1, X2] = BEM_StokesPartMeshLine(Xi, Xf, N)
    % INPUTS
    % 
    % Xi : 1st point of the line segment
    %
    % Xf : last point of the line segment
    % 
    % N : Number of elements along the line

    %{
        x coordinate of the initial point

    i : initial
    %}
    xi = Xi(1);

    % y coordinate of the initial point
    yi = Xi(2);

    %{
        x coordinate of the final point

    f : final
    %}
    xf = Xf(1);

    % y coordinate of the final point
    yf = Xf(2);

    % Initializing the matrix storing all the points on the line segment
    X_All = zeros(N + 1, 2);

    % x coordinates of the points
    X_All(:, 1) = linspace(xi, xf, N + 1);

    % x coordinates of the points
    X_All(:, 2) = linspace(yi, yf, N + 1);

    % 1st points of the elements
    X1 = X_All(1 : end - 1, :);

    % 2nd points of the elements
    X2 = X_All(2 : end, :);
end