function [result] = forward_kolmogorov(model,f,epsilon)
applyBoundaryCondition(model,"dirichlet","edge",[1,12],"u",0);
ccoeff = @(location, state) [0 -((location.y).^2)/2 ((location.y).^2)/2 epsilon];
specifyCoefficients(model,"m",0,"d",0,"c",ccoeff,"a",0,"f",f);
generateMesh(model);
result = solvepde(model); 
end