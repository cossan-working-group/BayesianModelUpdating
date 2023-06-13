function [proprnd] = proposal_rnd(CurrentSample, Tuning_mcmc, NumberOfChains, BoxFunction)

% This is a modified Proposal random number generator from the Normal
% Proposal distribution. The normrnd function is now multiplied by the Box
% function such that if the generated proposal sample values fall outide
% the range of values as stipulated by the Prior distribution, the function
% returns a 0 immediately and the sample value is rejected. This serves to
% prevent any values which fall outside the range of Prior values from 
% being accepted. 

proprnd_nominal = normrnd(CurrentSample,Tuning_mcmc,NumberOfChains,1);

% Initiate the storing of the array of Proposal sample values with an empty array:
proprnd = zeros(NumberOfChains,1);

% To store each array value of proprnd_nominal:
for i = 1:NumberOfChains
proprnd(i) = proprnd_nominal(i) .* BoxFunction(proprnd_nominal(i)); 
end

end