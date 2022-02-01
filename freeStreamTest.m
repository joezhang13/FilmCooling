function freeStreamTest(mesh)
%This function implements the free-stream preservation test. The initial
%condition and all the boundary conditions are set to the free-stream
%state. The residuals in the first 500 time steps are checked. 
%INPUT:
%  mesh = mesh structure

%Constants
CFL = 0.9; S0 = 8;
gamma = 1.4; RGas = 1; 
MInf = 0.25; rhoInf = 1; uInf = 1; fInf = 1; 
cInf = uInf / MInf; pInf = rhoInf * cInf^2 / gamma;
EInf = pInf / (rhoInf * (gamma - 1)) + 0.5 * uInf^2;
vPlenum = 0.05; etaTemp = 0.8; etaPres = 1.1;
TInf = pInf / (rhoInf * RGas); TPlenum = etaTemp * TInf;
pPlenum = etaPres * pInf; rhoPlenum = pPlenum / (RGas * TPlenum);
EPlenum = pPlenum / (rhoPlenum * (gamma - 1)) + 0.5 * vPlenum^2;
TOL = 10^(-5); 
Nt = 500;  %maximum number of iterations

%Read the mesh and set the initial condition
N = mesh.nNode;
nodes = mesh.Node;
Ne = mesh.nElem;
Elem = mesh.Elem;
B = mesh.B;
[IE,BE]=edgehash(Elem);
Ni = length(IE);
Nb = length(BE);
u = zeros(Ne, 5);
u(:, 1) = rhoInf; u(:, 2) = rhoInf * uInf;
u(:, 4) = rhoInf * EInf; u(:, 5) = rhoInf * fInf;

for n = 1 : Nt
    
    %Initialize the residual R_i and wave speed s_tot,i
    R = zeros(Ne, 5);
    stot = zeros(Ne, 1);

    %Loop over interior faces
    for i = 1 : Ni
        %Identify 2 cells (L, R) adjacent to face i
        neL = IE(i, 3);
        neR = IE(i, 4);
        uL = u(neL, :);
        uR = u(neR, :);

        %Determine the unit normal vector (L to R) and the length of edge
        elemL = Elem(neL, :);
        n1 = IE(i, 1);
        n2 = IE(i, 2);
        [nvec, deltaL] = normal(elemL, n1, n2, nodes);

        %Compute the flux and the wave speed using the flux function
        %and update the R and s_tot in each cell
        [F,smag] = FluxFunction(uL,uR,nvec);
        R(neL, :) = R(neL, :) + F' * deltaL; 
        R(neR, :) = R(neR, :) - F' * deltaL; 
        stot(neL) = stot(neL) + smag * deltaL;
        stot(neR) = stot(neR) + smag * deltaL;
    end

    %Loop over boundary faces
    for j = 1 : B.nbfgrp
        %Free-stream test
        for i = 1 : B.nbface(j)
            n1 = B.nodes{j}(i, 1);
            n2 = B.nodes{j}(i, 2);
            % ne = find((BE(:,1)==n1 & BE(:,2)==n2) | (BE(:,1)==n2 & BE(:,2)==n1));
            neIndex = (BE(:,1)==n1 & BE(:,2)==n2) | (BE(:,1)==n2 & BE(:,2)==n1);
            ne = BE(neIndex, 3);
            elemB = Elem(ne, :);
            uP = u(ne, :);
            %apply free-stream BC
            uB = rhoInf * [1, uInf, 0, EInf, fInf];
            [nvec, deltaL] = normal(elemB, n1, n2, nodes);
            [F, smag] = FluxFunction(uP, uB, nvec);
            R(ne, :) = R(ne, :) + F' * deltaL;
            stot(ne) = stot(ne) + smag * deltaL;
        end
    end

    %Add source term
    for i = 1 : Ne
        elem = Elem(i, :);
        ue = u(i, :);
        S = source(S0, ue);
        A = cellArea(elem, nodes);
        R(i, :) = R(i, :) + A * S;
    end

    %Compute the L1 norm of residuals
    Res = reshape(R', 5 * Ne, 1);
    L1(n) = norm(Res, 1);
    disp(['Iteration ', num2str(n),'. L1 norm of the residual: ', num2str(L1(n))]);

    %Update the state
    for i = 1 : Ne
        dt_A = 2 * CFL / stot(i);
        u(i, :) = u(i, :) - dt_A * R(i, :);
    end
    
end

%Plot the L1 norm of residuals at each time step
figure;
plot(1 : Nt, L1);
xlabel('number of time steps');
ylabel('L1 norm of residuals');
title('Free-stream preservation test');

end