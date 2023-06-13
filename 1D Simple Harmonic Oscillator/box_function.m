function [BoxFunction] = box_function(theta, parameters)

% The Box Function serves as an indicator function such that should the 
% Candidate samples fall outside the range of values of the Uniform Prior, 
% the function returns a 0 and returns a 1 if the values of the Candidate
% samples fall within the range of values of the Uniform Prior.

   
if theta < parameters(1) || theta > parameters(2)

BoxFunction = 0;
    
else
    
BoxFunction = 1;
    
end
    
end