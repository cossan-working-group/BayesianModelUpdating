function [samples_fT_D, log_fD] = TMCMCsampler(varargin)
%% Transitional Markov Chain Monte Carlo sampler
%
% This program implements a method described in:
% Ching, J. and Chen, Y. (2007). "Transitional Markov Chain Monte Carlo
% Method for Bayesian Model Updating, Model Class Selection, and Model
% Averaging." J. Eng. Mech., 133(7), 816-832.
%
% Usage:
% [samples_fT_D, fD] = tmcmc_v1(fD_T, fT, sample_from_fT, N);
%
% where:
%
% inputs:
% log_fD_T       = function handle of log(fD_T(t)), Loglikelihood
%
% fT             = function handle of fT(t), Prior PDF
%
% sample_from_fT = handle to a function that samples from of fT(t),
% Sampling rule function from Prior PDF
%
% nsamples              = number of samples of fT_D, Posterior, to generate
%
% outputs:
% samples_fT_D   = samples of fT_D (N x D) = samples from Posterior
% distribution
%
% log_fD         = log(evidence) = log(normalization constant)

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
   
%% Obtain N samples from the prior pdf f(T)
j      = 0;                   % Initialise loop for the transitional likelihood
thetaj = prior_rnd(nsamples);   % theta0 = N x D
pj     = 0;                   % p0 = 0 (initial tempering parameter)
Dimensions = size(thetaj, 2); % size of the vector theta

%% Initialization of matrices and vectors
thetaj1   = zeros(nsamples, Dimensions);
%log_fD_T_thetaj = zeros(nsamples,1);

%% Main loop
while pj < 1    
    j = j+1;
    
    %% Calculate the tempering parameter p(j+1):
    for l = 1:nsamples
        log_fD_T_thetaj(l) = loglikelihood(thetaj(l,:));
    end
    if any(isinf(log_fD_T_thetaj))
        error('The prior distribution is too far from the true region');
    end
    pj1 = calculate_pj1(log_fD_T_thetaj, pj);
    fprintf('TMCMC: Iteration j = %2d, pj1 = %f\n', j, pj1);
    
    %% Compute the plausibility weight for each sample wrt f_{j+1}
    fprintf('Computing the weights ...\n');
    % wj     = fD_T(thetaj).^(pj1-pj);         % N x 1 (eq 12)
    a       = (pj1-pj)*log_fD_T_thetaj;
    wj      = exp(a);
    wj_norm = wj./sum(wj);                % normalization of the weights
    
    %% Compute S(j) = E[w{j}] (eq 15)
    S(j) = mean(wj);
    
    %% Do the resampling step to obtain N samples from f_{j+1}(theta) and
    % then perform Metropolis-Hastings on each of these samples using as a
    % stationary PDF "fj1"
    % fj1 = @(t) fT(t).*log_fD_T(t).^pj1;   % stationary PDF (eq 11) f_{j+1}(theta)
    log_posterior = @(t) log(priorpdf(t)) + pj1*loglikelihood(t);
    

    % and using as proposal PDF a Gaussian centered at thetaj(idx,:) and
    % with covariance matrix equal to an scaled version of the covariance
    % matrix of fj1:
    
    % weighted mean
    mu = zeros(1, Dimensions);
    for l = 1:nsamples
        mu = mu + wj_norm(l)*thetaj(l,:); % 1 x N
    end
    
    % scaled covariance matrix of fj1 (eq 17)
    cov_gauss = zeros(Dimensions);
    for k = 1:nsamples
        % this formula is slightly different to eq 17 (the transpose)
        % because of the size of the vectors)m and because Ching and Chen
        % forgot to normalize the weight wj:
        tk_mu = thetaj(k,:) - mu;
        cov_gauss = cov_gauss + wj_norm(k)*(tk_mu'*tk_mu);
    end
    cov_gauss = beta^2 * cov_gauss;
    assert(~isinf(cond(cov_gauss)),'Something is wrong with the likelihood.')
    
    % Define the Proposal distribution:
    proppdf = @(x,y) prop_pdf(x, y, cov_gauss, priorpdf); %q(x,y) = q(x|y).
    proprnd = @(x)   prop_rnd(x, cov_gauss, priorpdf);   
    
    %% During the last iteration we require to do a better burnin in order
    % to guarantee the quality of the samples:
    if pj1 == 1
        burnin = lastBurnin;
    end;
    
    %% Start N different Markov chains
    fprintf('Markov chains ...\n\n');
    idx = randsample(nsamples, nsamples, true, wj_norm);
    for i = 1:nsamples      % For parallel, type: parfor
        %% Sample one point with probability wj_norm
        
        % smpl = mhsample(start, nsamples,
        %                'pdf', pdf, 'proppdf', proppdf, 'proprnd', proprnd);
        % start = row vector containing the start value of the Markov Chain,
        % nsamples = number of samples to be generated
        [thetaj1(i,:), acceptance_rate] = mhsample(thetaj(idx(i), :), 1, ...
            'logpdf',  log_posterior, ...
            'proppdf', proppdf, ...
            'proprnd', proprnd, ...
            'thin',    3,       ...
            'burnin',  burnin);
        % According to Cheung and Beck (2009) - Bayesian model updating ...,
        % the initial samples from reweighting and the resample of samples of
        % fj, in general, do not exactly follow fj1, so that the Markov
        % chains must "burn-in" before samples follow fj1, requiring a large
        % amount of samples to be generated for each level.
        
        %% Adjust the acceptance rate (optimal = 23%)
        % See: http://www.dms.umontreal.ca/~bedard/Beyond_234.pdf
        %{
      if acceptance_rate < 0.3
         % Many rejections means an inefficient chain (wasted computation
         %time), decrease the variance
         beta = 0.99*beta;
      elseif acceptance_rate > 0.5
         % High acceptance rate: Proposed jumps are very close to current
         % location, increase the variance
         beta = 1.01*beta;
      end
        %}
    end
    fprintf('\n');
    
    %% Prepare for the next iteration
    thetaj = thetaj1;
    pj     = pj1;
end

% TMCMC provides N samples distributed according to the Posterior distribution, f(T|D)
samples_fT_D = thetaj;

% estimation of f(D) -- this is the normalization constant in Bayes
log_fD = sum(log(S(1:j)));

return; % End


%% Calculate the tempering parameter p(j+1)
function pj1 = calculate_pj1(log_fD_T_thetaj, pj)
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

% wj = @(e) fD_T_thetaj^e; % N x 1
% Note the following trick in order to calculate e:
% Take into account that e>=0
wj = @(e) exp(abs(e)*log_fD_T_thetaj); % N x 1
%fmin = @(e) std(wj(e))/mean(wj(e)) - threshold;
fmin = @(e) std(wj(e)) - threshold*mean(wj(e)) + realmin;
e = abs(fzero(fmin, 0)); % e is >= 0, and fmin is an even function
if isnan(e)
    error('There is an error finding e');
end

pj1 = min(1, pj + e);

return; % End

function proppdf = prop_pdf(x, mu, covmat, box)
% This is the Proposal PDF for the Markov Chain.

% Box function is the Prior PDF in the feasible region. 
% So if a point is out of bounds, this function will
% return 0.

proppdf = mvnpdf(x, mu, covmat).*box(x); %q(x,y) = q(x|y).

return;


function proprnd = prop_rnd(mu, covmat, box)
% Sampling from the proposal PDF for the Markov Chain.

while true
    proprnd = mvnrnd(mu, covmat, 1);
    if box(proprnd)
        % The box function is the Prior PDF in the feasible region.
        % If a point is out of bounds, this function will return 0 = false.
        break;
    end
end

return
