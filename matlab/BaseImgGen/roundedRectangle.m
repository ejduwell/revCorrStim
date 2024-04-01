function outMat = roundedRectangle(x1, y1, x2, y2, r, np)
    % Define the number of points for the arcs
    %np = 100; % Increase for smoother curves

    % Define the angles for the four corners
    angles1 = linspace(pi, 3*pi/2, np);
    angles2 = linspace(3*pi/2, 2*pi, np);
    angles3 = linspace(0, pi/2, np);
    angles4 = linspace(pi/2, pi, np);

    % Calculate the arc points for each corner
    arc1_x = r*cos(angles1) + x1 + r;
    arc1_y = r*sin(angles1) + y1 + r;

    arc2_x = r*cos(angles2) + x2 - r;
    arc2_y = r*sin(angles2) + y1 + r;

    arc3_x = r*cos(angles3) + x2 - r;
    arc3_y = r*sin(angles3) + y2 - r;

    arc4_x = r*cos(angles4) + x1 + r;
    arc4_y = r*sin(angles4) + y2 - r;

    % Combine the points to form the outline
    x = [arc1_x, arc2_x, arc3_x, arc4_x, arc1_x(1)];
    y = [arc1_y, arc2_y, arc3_y, arc4_y, arc1_y(1)];

    outMat=horzcat(x',y');

    % Plot the rounded rectangle for visualization
%     plot(x, y, 'LineWidth', 2);
%     axis equal;
%     grid on;
%     xlabel('X');
%     ylabel('Y');
%     title('Rounded Rectangle');
end
