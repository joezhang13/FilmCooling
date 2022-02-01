function u = FVSolver(mesh, S0, CFL, u0)
%This function solves the governing equations (the Euler equations and the 
%fuel concentration equation) using first-order finete volume method
%INPUTS: 
%  mesh = mesh data structure
%  S0 = reaction rate amplitude
%  CFL = CFL number
%  u0 = initial state (free-stream state is used if not specified)
%OUTPUT:
%  u = state vectors for all elements (rho, rho*u, rho*v, rho*E, rho*f)

%Constants
gamma = 1.4; RGas = 1; 
MInf = 0.25; rhoInf = 1; uInf = 1; fInf = 1; 
cInf = uInf / MInf; pInf = rhoInf * cInf^2 / gamma;
EInf = pInf / (rhoInf * (gamma - 1)) + 0.5 * uInf^2;
vPlenum = 0.05; etaTemp = 0.8; etaPres = 1.1;
TInf = pInf / (rhoInf * RGas); TPlenum = etaTemp * TInf;
pPlenum = etaPres * pInf; rhoPlenum = pPlenum / (RGas * TPlenum);
EPlenum = pPlenum / (rhoPlenum * (gamma - 1)) + 0.5 * vPlenum^2;
TOL = 10^(-5); 
Nt = 25000;  %maximum number of iterations

%Read the mesh and set the initial condition
N = mesh.nNode;
nodes = mesh.Node;
Ne = mesh.nElem;
Elem = mesh.Elem;
B = mesh.B;
[IE,BE]=edgehash(Elem);
Ni = length(IE);
Nb = length(BE);
if nargin == 4
    u = u0;
elseif nargin == 3
    u = zeros(Ne, 5);
    u(:, 1) = rhoInf; u(:, 2) = rhoInf * uInf;
    u(:, 4) = rhoInf * EInf; u(:, 5) = rhoInf * fInf;
    disp('The default free-stream initial condition is used.');
else
    error('Wrong input!');
end

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
        group = B.title{j};
        
        %Free-stream boundaries
        if strcmp(group,'Inflow') || strcmp(group,'Top') || strcmp(group,'Outflow')
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

        %Plenum boundaries
        elseif strcmp(group,'Plenum')
            for i = 1 : B.nbface(j)
                n1 = B.nodes{j}(i, 1);
                n2 = B.nodes{j}(i, 2);
                % ne = find((BE(:,1)==n1 & BE(:,2)==n2) | (BE(:,1)==n2 & BE(:,2)==n1));
                neIndex = (BE(:,1)==n1 & BE(:,2)==n2) | (BE(:,1)==n2 & BE(:,2)==n1);
                ne = BE(neIndex, 3);
                elemB = Elem(ne, :);
                uP = u(ne, :);
                %apply plenum BC
                uB = rhoPlenum * [1, 0, vPlenum, EPlenum, 0];
                [nvec, deltaL] = normal(elemB, n1, n2, nodes);
                [F, smag] = FluxFunction(uP, uB, nvec);
                R(ne, :) = R(ne, :) + F' * deltaL;
                stot(ne) = stot(ne) + smag * deltaL;
            end

        %Wall boundaries
        elseif strcmp(group,'Wall1') || strcmp(group,'Wall2') || strcmp(group,'Wall3') || strcmp(group,'Wall4')
            for i = 1 : B.nbface(j)
                n1 = B.nodes{j}(i, 1);
                n2 = B.nodes{j}(i, 2);
                % ne = find((BE(:,1)==n1 & BE(:,2)==n2) | (BE(:,1)==n2 & BE(:,2)==n1));
                neIndex = (BE(:,1)==n1 & BE(:,2)==n2) | (BE(:,1)==n2 & BE(:,2)==n1);
                ne = BE(neIndex, 3);
                elemB = Elem(ne, :);
                uP = u(ne, :);
                %apply wall BC
                vPvec = [uP(2), uP(3)] / uP(1);
                [nvec, deltaL] = normal(elemB, n1, n2, nodes);
                un = dot(vPvec, nvec);
                vBvec = vPvec - un * nvec;
                rhoEP = uP(4); rhoP = uP(1);
                pB = (gamma - 1) * (rhoEP - 0.5 * rhoP * norm(vBvec)^2);
                F = [0; pB * nvec(1); pB * nvec(2); 0; 0];
                R(ne, :) = R(ne, :) + F' * deltaL;
                cB = (pB * gamma / rhoP)^0.5;
                l(1) = un + cB; l(2) = un - cB; l(3) = un; l = abs(l);
                smag = max(l);
                stot(ne) = stot(ne) + smag * deltaL;
            end
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

    %Check tolerance
    Res = reshape(R', 5 * Ne, 1);
    L1(n) = norm(Res, 1);
    disp(['Iteration ', num2str(n),'. L1 norm of the residual: ', num2str(L1(n))]);
    if L1(n) < TOL
        disp(['L1 residual norm drops below ', num2str(TOL), '. The running is ended!']);
        break;
    end

    %Update the state
    for i = 1 : Ne
        dt_A = 2 * CFL / stot(i);
        u(i, :) = u(i, :) - dt_A * R(i, :);
    end
    
    %Store the intermediate time steps for restart
    if ~mod(n, 1000)
        fileID = fopen('temp.txt', 'a');
        fprintf(fileID, '%18s\r\n', ['time step ', num2str(n)]);
        fprintf(fileID, '%12.8f %12.8f %12.8f %12.8f %12.8f\r\n', u');
        fclose(fileID);
    end
    
end

%Plot the L1 norm of residuals at each time step
figure;
plot(L1);
xlabel('number of time steps');
ylabel('L1 norm of residuals');

end