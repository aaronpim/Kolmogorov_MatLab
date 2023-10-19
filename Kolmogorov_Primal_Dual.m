epsilon = 1;
alpha = 0.1;
v = 1;
error_epsilon = 1e-8;
ud = @(location, state) (location.x<0.4).*(location.x>0.2).*(location.y>-0.3).*(location.y<-0.1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Domain Definition
% Define PDE model
model_forward = createpde(1);
model_adjoint = createpde(1);
% Define Geometry
gd = ...
[3, 3, 3, 3;...
 4, 4, 4, 4;...
 0, 0,-1,-1;...
 1, 1, 0, 0;...
 1, 1, 0, 0;...
 0, 0,-1,-1;...
 0,-v, 0,-v;...
 0,-v, 0,-v;...
 v, 0, v, 0;...
 v, 0, v, 0];
ns = [82, 82, 82, 82; 49, 50, 51, 52];
sf = 'R1+R2+R3+R4';
dl = decsg(gd,sf,ns);
% Assign Geometry to model
geometryFromEdges(model_forward,dl);
geometryFromEdges(model_adjoint,dl);
% Assign zero velocity boundary conditions
applyBoundaryCondition(model_forward,"dirichlet","edge",[3,4,5,6],"u",0);
applyBoundaryCondition(model_adjoint,"dirichlet","edge",[3,4,5,6],"u",0);
% Assign zero inflow and outflow conditions
applyBoundaryCondition(model_forward,"dirichlet","edge",[1,12],"u",0);
applyBoundaryCondition(model_adjoint,"dirichlet","edge",[2,11],"u",0);
% Generate mesh
generateMesh(model_forward);
generateMesh(model_adjoint);
%% Initial Dual solve
% Assign coefficients
syms z(x,y)
coeffs_adjoint = pdeCoefficients(-y.*diff(z,x) -epsilon*diff(z,y,y),z);
specifyCoefficients(model_adjoint,"m",0,"d",0,"c",coeffs_adjoint.c,"a",0,"f",ud);
results_adjoint = solvepde(model_adjoint);
%% Test forward solve
% Interpolation function
Err = error_epsilon + 1;
count = 0;
syms u(x,y)
coeffs_forward = pdeCoefficients(y.*diff(u,x) -epsilon*diff(u,y,y),u);
Err_vec = [];
while Err > error_epsilon
    if count > 1
        u_old = results_forward.NodalSolution;
    end
    f = @(location, state) transpose(max(interpolateSolution(results_adjoint,location.x,location.y),0))*(1/alpha);
  
    specifyCoefficients(model_forward,"m",0,"d",0,"c",coeffs_forward.c,"a",0,"f",f);
    results_forward = solvepde(model_forward);

    u = @(location, state) ud(location, state) - transpose(interpolateSolution(results_adjoint,location.x,location.y));
    specifyCoefficients(model_adjoint,"m",0,"d",0,"c",coeffs_adjoint.c,"a",0,"f",u);
    results_adjoint = solvepde(model_adjoint);
    if count > 1
        du = (results_forward.NodalSolution - u_old).^2;
        X  = results_forward.Mesh.Nodes(1,:);
        V  = results_forward.Mesh.Nodes(2,:);
        F = scatteredInterpolant(X',V',du);
        Err = integral2(@(x,y) F(x,y), -1, 1, -v, v);
        Err_vec = [Err_vec,Err];
    end
    count = count + 1;
end