function A = cellArea(elem, nodes)
%This function calculates the area of an element
%INPUTS:
%  elem = index of the element
%  nodes = array of node coordinates
%OUTPUT:
%  A = area of the element

tvec1 = [nodes(elem(2), :) - nodes(elem(1), :), 0];
tvec2 = [nodes(elem(3), :) - nodes(elem(1), :), 0];
A = 0.5 * abs(tvec1(1) * tvec2(2) - tvec1(2) * tvec2(1));

end