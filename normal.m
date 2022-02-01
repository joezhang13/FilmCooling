function [nvec, deltaL] = normal(elem, n1, n2, nodes) 
%This function calculates the normal vector and the length of an edge on a 
%given element
%INPUTS:
%  elem = index of the element
%  n1 = index of the first node
%  n2 = index of the second node
%  nodes = array of node coordinates
%OUTPUTS:
%  nvec = the normal vector (nx, ny)
%  deltaL = the length of the edge

%Define the edge vector clockwise
if (n1 == elem(1) && n2 == elem(2)) ||... 
   (n1 == elem(2) && n2 == elem(3)) ||... 
   (n1 == elem(3) && n2 == elem(1))
    tvec = nodes(n1, :) - nodes(n2, :);
else
    tvec = nodes(n2, :) - nodes(n1, :);
end
tvec = tvec';  %column vector
deltaL = norm(tvec);

%Rotate the edge vector in counter-clockwise direction 
%to get the normal vector 
Rot = [0, -1; 1, 0];
nvec = Rot * tvec / deltaL;
nvec = nvec';  %row vector

% %Plot for debugging
% vert = zeros(3, 2);
% vert(1, :) = nodes(elem(1), :);
% vert(2, :) = nodes(elem(2), :);
% vert(3, :) = nodes(elem(3), :);
% plot(vert([1, 2], 1), vert([1, 2], 2), 'k');
% hold on
% plot(vert([2, 3], 1), vert([2, 3], 2), 'k');
% plot(vert([3, 1], 1), vert([3, 1], 2), 'k');
% mid = 0.5 * (nodes(n1, :) + nodes(n2, :));
% plot([mid(1), mid(1) + nvec(1) * deltaL],... 
%      [mid(2), mid(2) + nvec(2) * deltaL], 'r');

end