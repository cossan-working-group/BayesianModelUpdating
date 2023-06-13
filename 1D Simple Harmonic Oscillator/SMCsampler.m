classdef SMCsampler
    
    properties  % Defining the inputs/characteristics of the class:
        
    % The following fields are defined as a {double} class as they are
    % numeric data variables. 
        
    % Note: {double} is the default numeric data type (class) in MATLAB
     
        N{double}               % No. of samples at each iteration
        theta{double}           % Initial samples 
        prior                   % Prior function
        log_likelihood          % Log Likelihood function
        K{double}               % No. of sampler iterations
        D{double}               % No. of dimensions
        q_cov{double}           % Proposal covariance
        logw{double}            % Log weights
        mean_estimate{double}   % Mean estimates of theta
        var_estimate{double}    % Variance estimates of theta
    end
    
    
methods 
    
function obj = SMCsampler(varargin)

% Parse the information in the name/value pairs: 
pnames = {'nsamples','prior_values','prior','loglikelihood',...
      'no_iterations','prop_covariance'};

dflts =  {[],[],[],[],[],[]}; % define default values
        
[obj.N,obj.theta,obj.prior,obj.log_likelihood,obj.K,obj.q_cov] = ...
         internal.stats.parseArgs(pnames, dflts, varargin{:});
                         
obj.D = size(obj.theta,2); % This defines the no. of dimensions of theta.
             
%% Find weights of the Prior samples:
        
% The weight of the prior samples is simply defined by the Likelihood function [Eq. (18)]:

             logw = zeros(obj.N,1);
             
               for i = 1:obj.N
                   logw(i) = obj.log_likelihood(obj.theta(i,:));
               end
          
             obj.logw = logw;
           
           
end


%% Log of the Posterior (Target) PDF:
        
% The Log_posterior is defined as the sum of the Log_prior and the 
% Log_likelihood functions:
        
function log_posterior = log_posterior(obj, theta)
log_posterior = obj.log_likelihood(theta) + log(obj.prior(theta));
end

%% Normalise weights function according to [Eq. (25)]:
        
        function norm_w = normalise_weights(obj, logw)
            w = exp(logw);
            norm_w = w ./ sum(w); 
        end

%% Re-sampling function:
       
        function [thetaNew, wn_new] = resample(obj, theta, wn)
           
            i_new = randsample(obj.N,obj.N,true,wn);
            
            % This is to reset the weights to 1/N during the Re-sampling process:
            wn_new = ones(obj.N,1) ./ obj.N; 
            thetaNew = theta(i_new,:);
            
        end

%% Estimate some quantities of interest:
        
        function [mean, variance] = estimate(obj, theta, wn)
           
            % Estimate the mean:
            mean = zeros(obj.D,1);
            for i = 1:obj.D
                mean(i) = sum( theta(:,i) .* wn );
            end

            % Estimate the variance:
            theta = theta - mean';
            variance = zeros(obj.D,1);
            for i = 1:obj.D
                variance(i) = sum(( theta(:,i).^2 ) .* wn);
            end
        end

%% Proposal function to generate proposal samples:
        
        function theta_proposed = propose_sample(obj, theta)
            if obj.D == 1
                
                while true
                    theta_proposed = theta + sqrt((obj.q_cov)) * randn();
                    if obj.prior(theta_proposed)
                        % For any problem, box is the function in the feasible region. So if
                        % a point is out of bounds set by Prior, this function will return 0 = false
                        break;
                    end
                end
            else
                while true
                    theta_proposed = mvnrnd(theta, obj.q_cov, 1);
                    if obj.prior(theta_proposed)
                        % For any problem, box is the function in the feasible region. So if
                        % a point is out of bounds set by Prior, this function will return 0 = false
                        break;
                    end
                end
            end
        end
 
%% Run the SMC sampler:

        function obj = generate_samples(obj)
            
        % Initialise the sampler vairables:
        logw_1 = obj.logw;
        theta_1 = obj.theta;
        theta_new = zeros(obj.N,obj.D);
        logw_new = zeros(obj.N,1);
        mean_estimate_1 = zeros(obj.K,obj.D);
        var_estimate_1 = zeros(obj.K,obj.D);

        % Main sampling loop: 
        for k = 1:obj.K
            wn = obj.normalise_weights(logw_1);
            [mean_estimate_1(k,:), var_estimate_1(k,:)] = obj.estimate(theta_1,wn);

           
            Neff = 1/sum(wn.^2); % Calculate the effective sample size using [Eq. (22)]
            
            % Resample if effective sample size is below threshold N/2:
            if Neff < obj.N/2
                [theta_1,wn] = obj.resample(theta_1, wn);
                logw_1 = log(wn);
            end
        
            % Generate new proposed samples and find new weights of each propsoed samples using [Eq. (25)]:
            for i = 1:obj.N
                theta_new(i,:) = obj.propose_sample(theta_1(i,:));
                logw_new(i) = logw_1(i) + obj.log_posterior(theta_new(i,:)) - obj.log_posterior(theta_1(i,:));
            end
            theta_1 = theta_new;  % Update samples
            logw_1 = logw_new;    % Update log weights
        end
        
        % Report final samples, log weights and quantities being estimated:
         obj.theta = theta_1;
         obj.logw = logw_1;
         obj.mean_estimate = mean_estimate_1;
         obj.var_estimate = var_estimate_1;
    
        end
    end
end