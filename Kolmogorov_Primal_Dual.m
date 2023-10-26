warning('off','all')

%% Parameters
epsilon = 1;            % Rate of diffusion
alpha = 0.1;              % Weight of source term
v = 1;                  % Maximum velocity
error_epsilon = 1e-8;   % Maximum relative error
%% Target function
ud = @(location, state) (location.x<0.4).*(location.x>0.2);
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
hmax = model_forward.Mesh.MaxElementSize;
hmin = model_forward.Mesh.MinElementSize;

% Refine mesh
generateMesh(model_forward,"Hmax",hmax/3,"Hmin",hmin/3);
generateMesh(model_adjoint,"Hmax",hmax/3,"Hmin",hmin/3);

%% Initial Dual solve

% Assign Dual equation coefficients
syms z(x,y)
coeffs_adjoint = pdeCoefficients(-y.*diff(z,x) -epsilon*diff(z,y,y),z);
specifyCoefficients(model_adjoint,"m",0,"d",0,"c",coeffs_adjoint.c,"a",0,"f",ud);

% Solve adjoint equation with u=0
results_adjoint = solvepde(model_adjoint);

%% Initialise Error vector and count variable
Err = error_epsilon + 1;
Err_vec = [];
count = 0;

% Assign Primal equation coefficients
syms u(x,y)
coeffs_forward = pdeCoefficients(y.*diff(u,x) -epsilon*diff(u,y,y),u);

while Err > error_epsilon

    % Record the previous result
    if count > 1
        u_old = results_forward.NodalSolution;
    end

    % Take the projection of the solution to the adjoint equation in only
    % positive values. Then scale by 1/alpha
    f = @(location, state) transpose(max(interpolateSolution(results_adjoint,location.x,location.y),0))*(1/alpha);

    % Redefine source term in forward model using the above
    specifyCoefficients(model_forward,"m",0,"d",0,"c",coeffs_forward.c,"a",0,"f",f);

    % Solve forward equation
    results_forward = solvepde(model_forward);
    
    % Define difference between solution and target solution
    u = @(location, state) ud(location, state) - transpose(interpolateSolution(results_adjoint,location.x,location.y));

    % Redefine source term in dual model using the above
    specifyCoefficients(model_adjoint,"m",0,"d",0,"c",coeffs_adjoint.c,"a",0,"f",u);

    % Solve the dual eqation
    results_adjoint = solvepde(model_adjoint);

    %% Calculate relative L2 error between one iteration and the next.
    if count > 1
        du = ((results_forward.NodalSolution - u_old)./results_forward.NodalSolution).^2;
        du(isnan(du)) =0;
        X  = results_forward.Mesh.Nodes(1,:);
        V  = results_forward.Mesh.Nodes(2,:);
        F = scatteredInterpolant(X',V',du);
        Err = sqrt(integral2(@(x,y) F(x,y), -1, 1, -v, v));
        Err_vec = [Err_vec,Err];
    end
    count = count + 1;
end
%% Error Plot
semilogy(Err_vec);
xlabel('No. of Primal-Dual iterations')
ylabel('Relative Error: Particle Density')
%% Dose Plot
% Calculate spatial sample points
N = 100;
x_d = linspace(-1,1,N+1);
dose = zeros(N,1);

X = results_forward.Mesh.Nodes(1,:);
V = results_forward.Mesh.Nodes(2,:);
U = scatteredInterpolant(X',V',results_forward.NodalSolution);
for i = 1:N 
    dose(i) = integral2(@(x,y) U(x,y), x_d(i), x_d(i+1), -v, v);
end
figure
plot(x_d(1:end-1),dose./max(dose),'-k')
hold on
plot([0.2,0.2],[0,1],'--r'); plot([0.4,0.4],[0,1],'--r');
ylabel('Normalised Dose'); xlabel('Position')
