function [result] = forward_kolmogorov(model,f,epsilon)
applyBoundaryCondition(model,"dirichlet","edge",[1,12],"u",0);
syms u(x,y)
pdeeq = y.*diff(u,x) -epsilon*diff(u,y,y);
coeffs = pdeCoefficients(pdeeq,u);
specifyCoefficients(model,"m",0,"d",0,"c",coeffs.c,"a",0,"f",f);
generateMesh(model);
result = solvepde(model); 
end