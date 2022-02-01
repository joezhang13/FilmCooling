function mesh = newMesh(mesh0, alphaDegree)
%This function reads the original mesh struct, changes the jet angle
%alpha and output the new mesh struct. A plot of the new mesh will also
%be generated.
%INPUTS:
%  mesh0 = the original mesh with alpha equals to 45 degrees
%  alphaDegree = the jet angle alpha for the new mesh
%OUTPUT:
%  mesh = new mesh structure

mesh = mesh0;
N = mesh.nNode;
nodes = mesh.Node;
Elem=mesh.Elem;
[IE,BE]=edgehash(Elem);

Y0 = 0;
H = 0.5;
alpha = (alphaDegree / 180) * pi;

for i = 1 : N
    y = nodes(i, 2);
    %move the channal nodes
    if y < Y0 && y > -H
        delta = (Y0 - y) * (1 - cot(alpha));
        nodes(i, 1) = nodes(i, 1) + delta;
    %move the plenum nodes
    elseif y <= -H
        delta = H * (1 - cot(alpha));
        nodes(i, 1) = nodes(i, 1) + delta;
    end
end
mesh.Node = nodes;  %write to new mesh

%plot the new mesh
figure;
hold on
for i = 1 : length(IE)
    n1 = IE(i, 1);
    n2 = IE(i, 2);
    plot([nodes(n1, 1) nodes(n2, 1)], [nodes(n1, 2) nodes(n2, 2)], 'k');
end
for i = 1 : length(BE)
    n1 = BE(i, 1);
    n2 = BE(i, 2);
    plot([nodes(n1, 1) nodes(n2, 1)], [nodes(n1, 2) nodes(n2, 2)], 'b');
end
axis equal;
ylim([-1.5 1.5]);
xlim([0 6]);
title(['\alpha = ', num2str(alphaDegree), '\circ']);

end