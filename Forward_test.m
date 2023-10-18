[model] = domain_defn(1);
f = @(location, state) (location.x<0.4).*(location.x>0.2).*(location.y>-0.3).*(location.y<-0.1);
[result] = forward_kolmogorov(model,f,0.1);