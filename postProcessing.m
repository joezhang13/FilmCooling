function [TB, xB] = postProcessing(mesh, u)
%This function reads the converged results and generate the field plots of
%the Mach number, density, temperature and fuel conventration. The
%normalized temperature along the wall downstream of the cooling hole is
%returned.
%INPUTS:
%  mesh = mesh structure
%  u = state vectors for all elements (rho, rho*u, rho*v, rho*E, rho*f)
%OUTPUTS:
%  TB = normalized temperature T/T_Inf along the downstream wall
%  xB = x coordinates corresponding to TB along the downstream wall

%Constants
gamma = 1.4; RGas = 1;
MInf = 0.25; rhoInf = 1; uInf = 1;
cInf = uInf / MInf; pInf = rhoInf * cInf^2 / gamma;
TInf = pInf / (rhoInf * RGas);

%Read mesh and states
nodes = mesh.Node;
Elem = mesh.Elem;
nn = mesh.nNode;
ne = mesh.nElem;
B = mesh.B;
[~, BE] = edgehash(Elem);
rho = u(:, 1); uu = u(:, 2) ./ rho; vv = u(:, 3) ./ rho; 
E = u(:, 4) ./ rho; f = u(:, 5) ./ rho;

%Calculate the results
p = (gamma - 1) .* rho .* (E - 0.5 * (uu.^2 + vv.^2));
T = p ./ (RGas * rho);
c = (p * gamma ./ rho).^0.5; U = (uu.^2 + vv.^2).^0.5;
M = U ./ c;
for j = 1 : B.nbfgrp
    group = B.title{j};
    if strcmp(group,'Wall3')
        for i = 1 : B.nbface(j)
            n1 = B.nodes{j}(i, 1);
            n2 = B.nodes{j}(i, 2);
            neIndex = (BE(:,1)==n1 & BE(:,2)==n2) | (BE(:,1)==n2 & BE(:,2)==n1);
            ne = BE(neIndex, 3);
            elemB = Elem(ne, :);
            uP = u(ne, :);
            %apply wall BC
            vPvec = [uP(2), uP(3)] / uP(1);
            [nvec, ~] = normal(elemB, n1, n2, nodes);
            un = dot(vPvec, nvec);
            vBvec = vPvec - un * nvec;
            rhoEP = uP(4); rhoP = uP(1);
            pB = (gamma - 1) * (rhoEP - 0.5 * rhoP * norm(vBvec)^2);
            TB(i) = pB / (RGas * rhoP);
            xB(i) = 0.5 * (nodes(n1, 1) + nodes(n2, 1));
        end
    end
end
TB = TB / TInf;

%Field plots
fieldResults = [M, rho, T, f];
fieldTitles = {'Mach number', 'Density', 'Temperature', 'Fuel concentration'};
for i = 1 : 4
    val = fieldResults(:, i);
    figure;
    trisurf(Elem, nodes(:, 1), nodes(:, 2), zeros(nn, 1), val);
    view(2);
    axis equal
    colorbar;
    ylim([-1.5 1.5]);
    xlim([0 6]);
    title(fieldTitles{i});
    set(gca, 'position', [0.1, 0.05, 0.75, 0.95]);
end

end