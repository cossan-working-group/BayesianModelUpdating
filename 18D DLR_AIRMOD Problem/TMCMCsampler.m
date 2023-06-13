function [samples, log_evidence] = TMCMCsampler(varargin)
%% Transitional Markov Chain Monte Carlo sampler
%
% This program implements a method described in:
% Ching, J. and Chen, Y. (2007). "Transitional Markov Chain Monte Carlo
% Method for Bayesian Model Updating, Model Class Selection, and Model
% Averaging." J. Eng. Mech., 133(7), 816-832.
%
%
% inputs:
% loglikelihood  = Loglikelihood function handle
%
% priorpdf       = Prior PDF function handle;
%
% prior_rnd      = Random number generator from Prior function;
%
% nsamples       = number of samples to generate from Posterior;
%
% outputs:
% samples        = samples from Posterior distribution;
%
% log_evidence   = log(evidence) = log(normalization constant);

% ------------------------------------------------------------------------
% who                    when         observations
%--------------------------------------------------------------------------
% Diego Andres Alvarez   Jul-24-2013  First algorithm
%--------------------------------------------------------------------------
% Diego Andres Alvarez - daalvarez@unal.edu.co
% Edoardo Patelli      - edoardo.patelli@strath.ac.uk

% parse the information in the name/value pairs: 
pnames = {'nsamples','loglikelihood','priorpdf','priorrnd','burnin','lastburnin','beta'};

dflts =  {[],[],[],[],[],0,0.2}; % define default values
      
[nsamples,loglikelihood,priorpdf,prior_rnd,burnin,lastBurnin,beta] = ...
       internal.stats.parseArgs(pnames, dflts, varargin{:});
   
%% Obtain N samples from the Prior pdf:
j      = 0;                   % Initialise loop for the transitional likelihood
thetaj = prior_rnd(nsamples); % Generate initial samples from the Prior
pj     = 0;                   % p0 = 0 (initial value of tempering parameter, beta_j)
Dimensions = size(thetaj, 2); % This defines the no. of dimensions of theta.

%% Initialization of matrices and vectors
thetaj1   = zeros(nsamples, Dimensions);

%% Main loop
while pj < 1    
    j = j+1;
    
    %% Calculate the tempering parameter beta_{j+1}:
    
    % To calculate and store Likelihood values for each theta:
    for i = 1:nsamples
        log_likelihood(i) = loglikelihood(thetaj(i,:));
    end
    if any(isinf(log_likelihood))
        error('The prior distribution is too far from the true region');
    end
    
    % To calculate the tempering parameter (beta_j) according to [Eq. (17)]:
    pj1 = calculate_pj1(log_likelihood, pj);
    fprintf('TMCMC: Iteration j = %2d, pj1 = %f\n', j, pj1);
    
    %% Compute the normalised weights for each sample:
    
    fprintf('Computing the weights ...\n');
    a       = (pj1-pj)*log_likelihood;
    wj      = exp(a);
    wj_norm = wj./sum(wj); % Normalization of the weights [Eq. (18)]
    
    %% Compute S(j): The mean of the weights
    S(j) = mean(wj);
    
    %% Define the Transitional distribution, P^{j+1}:
    log_posterior = @(t) log(priorpdf(t)) + pj1*loglikelihood(t);
    

    % Define the proposal PDF a Gaussian centered at current sample (thetaj(idx,:)) and
    % with covariance matrix equal to an scaled version of the covariance
    % matrix of Transitional distribution P^{j+1}:
    
    % Weighted mean of theta:
    mu = zeros(1, Dimensions);
    for i = 1:nsamples
        mu = mu + wj_norm(i)*thetaj(i,:); % 1 x N
    end
    
    % Scaled covariance matrix [Eq. (19)]
    cov_gauss = zeros(Dimensions);
    for k = 1:nsamples
        tk_mu = thetaj(k,:) - mu;
        cov_gauss = cov_gauss + wj_norm(k)*(tk_mu'*tk_mu);
    end
    cov_gauss = beta^2 * cov_gauss;
    assert(~isinf(cond(cov_gauss)),'Something is wrong with the likelihood.')
    
    % Define the Proposal distribution:
    proppdf = @(x,y) prop_pdf(x, y, cov_gauss, priorpdf); %q(x,y) = q(x|y).
    proprnd = @(x)   prop_rnd(x, cov_gauss, priorpdf);   
    
    %% To set the burn-in for the last iteration (if needed)
    if pj1 == 1
        burnin = lastBurnin;
    end;
    
    %% Start N different Markov chains
    fprintf('Markov chains ...\n\n');
    idx = randsample(nsamples, nsamples, true, wj_norm);
    for i = 1:nsamples      % For parallel, type: parfor
        %% Sample one point via MH sampling with probability wj_norm

        [thetaj1(i,:), acceptance_rate] = mhsample(thetaj(idx(i), :), 1, ...
            'logpdf',  log_posterior, ...
            'proppdf', proppdf, ...
            'proprnd', proprnd, ...
            'thin',    3,       ...
            'burnin',  burnin);

    end
    fprintf('\n');
    
    %% Prepare for the next iteration
    thetaj = thetaj1;
    pj     = pj1;
end

% TMCMC provides N samples distributed according to the Posterior
% distribution:
samples = thetaj;

% Estimation of the Log(evidence):
log_evidence = sum(log(S(1:j)));

return; % End


%% Calculate the tempering parameter beta_{j+1}
function pj1 = calculate_pj1(log_likelihood, pj)
% find pj1 such that COV <= threshold, that is
%
%  std(wj)
% --------- <= threshold
%  mean(wj)
%
% here
% size(thetaj) = N x D,
% wj = fD_T(thetaj).^(pj1 - pj)
% e = pj1 - pj

threshold = 1; % 100% = threshold on the COV

% Note the following trick in order to calculate e:
% Take into account that e >= 0
wj = @(e) exp(abs(e)*log_likelihood); % N x 1
fmin = @(e) std(wj(e)) - threshold*mean(wj(e)) + realmin;
e = abs(fzero(fmin, 0)); % e is >= 0, and fmin is an even function
if isnan(e)
    error('There is an error finding e');
end

pj1 = min(1, pj + e);

return; % End

function proppdf = prop_pdf(candidate_sample, current_sample, covmat, box)
% This is the Proposal PDF for the Markov Chain.

% Box function is the Prior PDF in the feasible region. 
% So if a point is out of bounds, this function will
% return 0.

proppdf = mvnpdf(candidate_sample, current_sample, covmat).*box(candidate_sample); %q(x,y) = q(x|y).

return;


function proprnd = prop_rnd(current_sample, covmat, box)
% Sampling from the proposal PDF for the Markov Chain.

while true
    proprnd = mvnrnd(current_sample, covmat, 1);
    if box(proprnd)
        % The box function is the Prior PDF in the feasible region.
        % If a point is out of bounds, this function will return 0 = false.
        break;
    end
end

return
