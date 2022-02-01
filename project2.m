clear;

%Mesh
mesh45 = read_gri('mesh.gri');
mesh60 = newMesh(mesh45, 60);
mesh90 = newMesh(mesh45, 90);
mesh120 = newMesh(mesh45, 120);

%Free-stream preservation test
freeStreamTest(mesh45);

%Solve the problem with no reaction
CFL = 0.9; S0 = 0;
u45 = FVSolver(mesh45, S0, CFL);
u60 = FVSolver(mesh60, S0, CFL, u45);
u90 = FVSolver(mesh90, S0, CFL, u60);
u120 = FVSolver(mesh120, S0, CFL, u90);

%Solve the problem with reaction
S0 = 8;
u45S8 = FVSolver(mesh45, S0, CFL, u45);
u60S8 = FVSolver(mesh60, S0, CFL, u60);
u90S8 = FVSolver(mesh90, S0, CFL, u90);
u120S8 = FVSolver(mesh120, S0, CFL, u120);

%Post-processing

%After running the cases above, the saved results can be used for
%post-processing:

% load results.mat;

[TB45, xB45] = postProcessing(mesh45, u45);
[TB60, xB60] = postProcessing(mesh60, u60);
[TB90, xB90] = postProcessing(mesh90, u90);
[TB120, xB120] = postProcessing(mesh120, u120);
[TB45S8, xB45S8] = postProcessing(mesh45, u45S8);
[TB60S8, xB60S8] = postProcessing(mesh60, u60S8);
[TB90S8, xB90S8] = postProcessing(mesh90, u90S8);
[TB120S8, xB120S8] = postProcessing(mesh120, u120S8);
%normalized temperature for cases with no reaction
figure;
plot(xB45, TB45, 'LineWidth',2);
hold on
plot(xB60, TB60, 'LineWidth',2);
plot(xB90, TB90, 'LineWidth',2);
plot(xB120, TB120, 'LineWidth',2);
legend('\alpha = 45^\circ', '\alpha = 60^\circ', '\alpha = 90^\circ', '\alpha = 120^\circ');
xlabel('x');
ylabel('T/T_\infty');
%normalized temperature for cases with reaction
figure;
plot(xB45S8, TB45S8, 'LineWidth',2);
hold on
plot(xB60S8, TB60S8, 'LineWidth',2);
plot(xB90S8, TB90S8, 'LineWidth',2);
plot(xB120S8, TB120S8, 'LineWidth',2);
legend('\alpha = 45^\circ', '\alpha = 60^\circ', '\alpha = 90^\circ', '\alpha = 120^\circ');
xlabel('x');
ylabel('T/T_\infty');