function logl = log_likelihood(Thetas, Eigenvalues, standard_deviations, ModelHandle)
% Calculation of the log_likelihood for 2D Inverse Eigenvalue Problem:
%
% USAGE:
% logl = log_likelihood(Thetas, Eigenvalues, standard_deviations, ModelHandle)
%
% INPUTS:
% Thetas = Vector of epistemic parameters theta_1 and theta_2  [Nsamples x 2]
% Eigenvalues = Observations of Eigenvalues from the 2 x 2 Square Matrix  [Nobservations x 2]
% standard_deviations = Vector of the standard deviations of the 2D log likelihood function  [2 x 1]
% ModelHandle = the function handle of the model (see file "model.m")
%
% OUTPUTS:
% logl = loglikelihood function for the set of estimated values of theta_1
% and theta_2 as well as the Eigenvalues


%% Evaluate the model:
nchains=size(Thetas,1);
logl=zeros(nchains,1);

for n=1:nchains
modelOutput = ModelHandle(Thetas(n,:));

% Note: Details to the model can be found in the file: "model.m"

%% Compute the overall 2D log-likelihood:
    
logl_1 = - 0.5 * (1/standard_deviations(1))^2 *(Eigenvalues(:,1) - modelOutput(1))' * (Eigenvalues(:,1) - modelOutput(1));
    
% Note: logl_1 is the loglikelihood function in the theta_1 dimension.

logl_2 = - 0.5 * (1/standard_deviations(2))^2 *(Eigenvalues(:,2) - modelOutput(2))' * (Eigenvalues(:,2) - modelOutput(2));

% Note: logl_2 is the loglikelihood function in the theta_2 dimension.

logl(n) = logl_1 + logl_2;

% logl is the overall 2D loglikelihood function.
    
end