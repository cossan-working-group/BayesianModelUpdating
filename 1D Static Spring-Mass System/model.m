function [force] = model(stiffness,displacement)

%MODEL: 1-Degree of freedom mass-spring system 

force = - stiffness.*displacement;

end

