function logL = Airmod_log_p_D_theta(D, theta, Xnn_no_data, selected_modes)
% Calculation of the log_likelihood for the Airmod structure
%
% USAGE:
% logL = problemA_log_p_D_theta(D, theta)
%
% INPUTS:
% D     = experimental observations   nobs x dim_x
% theta = epistemic parameters        npar x dim_theta
%
% OUTPUTS:
% logL(i)  = loglikelihood for the set of parameters theta(i,:) and the
%            data D, i = 1, 2, ...npar.        logL = npar x 1

%%
npar = size(theta,1);  % number of thetas to evaluate
logL = zeros(npar,1);

% call the model once for all the samples

% load('metamodel_light','Xnn_no_data');

%% Set input values
Tinput = cell2struct(num2cell(theta(:,1:18)),Xnn_no_data(1).Cinputnames,2);
%% Map these simulations through the system
output = zeros(length(Tinput),length(selected_modes));
for imode = 1:size(D,2)
    Xout = Xnn_no_data(selected_modes(imode)).apply(Tinput);
    output(:,imode)  = Xout.getValues('CSnames',Xnn_no_data(selected_modes(imode)).Coutputnames(1));
end

%% compute loglikelihood of each sample
for ipar = 1:npar
    logL(ipar) = sum((p_x_theta_pdf(D, output(ipar,:), theta(ipar,:))));
    if isinf(logL(ipar))
        logL(ipar) = -1e10;
    end
end

end

%%
function p = p_x_theta_pdf(w, w_model, theta_i)
% x       = set of observations         nobs x dim_x
% theta_i = point in epistemic space    nsamples x dim_theta

fixed_variance = 0;

if ~fixed_variance
    epsilon_r = theta_i(:,19:end); % std multiplied by 1000000 for numerical reasons
else
    epsilon_r = 1e-4*ones(size(theta_i));
end
p = zeros(size(theta_i,1),size(w_model,2));
%% Estimate the PDF p_x_theta_pdf(x | theta)
% compute likelihood accordind to the GOCE paper
% p(i) = p_x_theta_pdf(x(i,:) | theta)
for i = 1:length(w_model)
    p(i) = (-1/2*nansum(((1-w(:,i).^2/w_model(i)^2)/epsilon_r(i)).^2,1));
end

end