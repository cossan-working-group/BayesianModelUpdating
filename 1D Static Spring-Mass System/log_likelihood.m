function logl = log_likelihood(stiffness, modelInput, measurements, standard_deviation, ModelHandle)
% Calculation of the log_likelihood for 1D Linear Static Problem:
%
% USAGE:
% logl = log_likelihood(stiffness, displacement, measurements, standard_deviation, ModelHandle)
%
% INPUTS:
% stiffness = epistemic parameter k          [Nsamples x 1]
% displacement = model input                 [Nobservations x 1]
% measurements = experimental observations   [Nobservations x 1]
% standard_deviation = the standard deviation of the log likelihood function  [scalar]
% ModelHandle = the function handle of the model (see file "model.m")
%
% OUTPUTS:
% logl = loglikelihood function for the set of estimated stiffness values k and
% the measurements 

    
%% Evaluate the model:
nchains=size(stiffness,1);
%ndims=size(stiffness,2);
logl=zeros(nchains,1);

for n=1:nchains
modelOutput = ModelHandle(stiffness(n,:),modelInput); 
  
% Note: Details to the model can be found in the file: "model.m"

%% Compute the log-likelihood:

logl(n) = -0.5 * (1/standard_deviation)^2 *(measurements - modelOutput)' * (measurements - modelOutput);
    
end

