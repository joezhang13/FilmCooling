function [F,smag] = FluxFunction(UL,UR,n)
    gamma = 1.4;
    gmi = gamma - 1.0;
    %Computes the Roe Flux in 2D
    %INPUTS:
    %   UL = L state vector (5 numbers = rho, rho*u, rho*v, rho*E, rho*f)
    %   UR = R state vector (5 numbers = rho, rho*u, rho*v, rho*E, rho*f)
    %   n  = unit spatial normal (2 numbers = nx, ny)
    %OUTPUTS:
    %   F = normal flux (5 numbers)
    %   smag = maximum wave speed magnitude
    
    %process left state
    rL = UL(1);
    uL = UL(2)/rL;
    vL = UL(3)/rL;
    unL = uL*n(1) + vL*n(2);
    qL = (UL(2)^2 + UL(3)^2)^(1/2) / rL;
    pL = (gamma-1)*(UL(4) - 0.5*rL*qL^2);
    if pL<0 || rL<0
        disp('Non-physical L state!');
        disp(UL);
    end
    rHL = UL(4) + pL;
    HL = rHL/rL;
    cL = (gamma*pL/rL)^0.5;
    fL = UL(5)/rL;      %fuel concentration for the L cell
    
    %left flux
    FL = zeros(4,1);
    FL(1) = rL*unL;
    FL(2) = UL(2)*unL + pL*n(1);
    FL(3) = UL(3)*unL + pL*n(2);
    FL(4) = rHL*unL;
    FL(5) = UL(5)*unL;      %fuel flux for the L cell
    
    %process right state
    rR = UR(1);
    uR = UR(2)/rR;
    vR = UR(3)/rR;
    unR = uR*n(1) + vR*n(2);
    qR = (UR(2)^2 + UR(3)^2)^(1/2) / rR;
    pR = (gamma-1)*(UR(4) - 0.5*rR*qR^2);
    if pR<0 || rR<0
        disp('Non-physical R state!');
        disp(UR);
    end
    rHR = UR(4) + pR;
    HR = rHR/rR;
    cR = (gamma*pR/rR)^0.5;
    fR = UR(5)/rR;      %fuel concentration for the R cell
    
    %right flux
    FR = zeros(4,1);
    FR(1) = rR*unR;
    FR(2) = UR(2)*unR + pR*n(1);
    FR(3) = UR(3)*unR + pR*n(2);
    FR(4) = rHR*unR;
    FR(5) = UR(5)*unR;      %fuel flux for the R cell
    
    %difference in states
    du = UR -UL;
    
    %Roe average;
    di = (rR/rL)^0.5;
    d1 = 1/(1+di);
    
    ui = (di*uR + uL)*d1;
    vi = (di*vR + vL)*d1;
    Hi = (di*HR + HL)*d1;
    
    af = 0.5*(ui^2 + vi^2);
    ucp = ui*n(1) + vi*n(2);
    c2 = gmi*(Hi - af);
    if c2<0
       disp('Non-physical Roe-average sound speed state!');
       c2 = -c2;
    end
    ci = c2^0.5;
    ci1 = 1.0/ci;
    
    %eigenvalues
    l = zeros(3,1);
    l(1) = ucp + ci;
    l(2) = ucp - ci;
    l(3) = ucp;
    
    %entropy fix
    epsilon = ci*0.02;
    for i = 1:3
        if abs(l(i)) < epsilon
            l(i) = 0.5*(epsilon + l(i)^2/epsilon);
        end
    end
    
    l = abs(l); l3 = l(3);
    
    %average and half-difference of 1st and 2nd eigs
    s1 = 0.5*(l(1) + l(2));
    s2 = 0.5*(l(1) - l(2));
    
    %left eigenvector product generators
    G1 = gmi*(af*du(1) - ui*du(2) - vi*du(3) + du(4));
    G2 = -ucp*du(1) + du(2)*n(1) + du(3)*n(2);
    
    %required functions of G1 and G2 (again, see Theory guide)
    C1 = G1*(s1 - l3)*ci1^2 + G2*s2*ci1;
    C2 = G1*s2*ci1 + G2*(s1 - l3);
    
    %flux assembly
    F = zeros(5,1);
    F(1) = 0.5*(FL(1) + FR(1)) - 0.5*(l3*du(1) + C1);
    F(2) = 0.5*(FL(2) + FR(2)) - 0.5*(l3*du(2) + C1*ui + C2*n(1));
    F(3) = 0.5*(FL(3) + FR(3)) - 0.5*(l3*du(3) + C1*vi + C2*n(2));
    F(4) = 0.5*(FL(4) + FR(4)) - 0.5*(l3*du(4) + C1*Hi + C2*ucp);
    F(5) = 0.5*(FL(5) + FR(5)) - 0.5*abs(ucp)*(rR*fR - rL*fL);  %numerical fuel flux calculated by a upwind method
    
    smag = max(l);
end