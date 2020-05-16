# Bayesian Model Updating: 

Bayesian Model Updating is a technique which casts the model updating problem in the form of a Bayesian Inference. There have been 3 popular advanced Monte Carlo sampling techniques which are adopted by researchers to address Bayesian Model Updating problems and make the necessary estimations of the epistemic parameter(s). These 3 techniques are:

* Markov Chain Monte Carlo [(MCMC)](https://doi.org/10.1093/biomet/57.1.97)
* Transitional Markov Chain Monte Carlo [(TMCMC)](https://doi.org/10.1061/(ASCE)0733-9399(2007)133:7(816))
* Sequential Monte Carlo [(SMC)](https://www.jstor.org/stable/3879283)

In this repository, 4 tutorials are presented to enable users to understand how the advanced Monte Carlo techniques are implemented in addressing various Bayesian Model Updating problems. The following tutorials are (in order of increasing difficulty):

* 1-Dimensional Linear Static System
* 1-Dimensional Simple Harmonic Oscillator
* 2-Dimensional Eigenvalue Problem
* 18-Dimensional DLR-AIRMOD Problem


## Tutorials:

### 1-Dimensional Static Spring-Mass System:

This tutorial presents a simple static Spring-Mass system. In this set-up, the spring is assumed to obey [Hooke's Law](http://latex.codecogs.com/svg.latex?F%3D-k%5Ccdot%7Bd%7D) model whereby the restoring force of the spring, F, is linearly proportional to the length of its displacement from rest length, d. The elasticity constant of the spring is k. This study, seeks to realize two objectives: 

1. To compare the estimation the epistemic parameter k between the samplers;

2. To compare the model updating results obtained through the use of MCMC, TMCMC, and SMC.

### 1-Dimensional Simple Harmonic Oscillator:

This tutorial presents a simple harmonic oscillator system. In this set-up, the natural oscillating frequency of the ocillator, F, obeys the [Simple Harmonic Frequency](http://latex.codecogs.com/svg.latex?F%3D%5Csqrt%7B%5Cfrac%7Bk%7D%7Bm%7D%7D) model whereby F is defined as the square-root of the ratio between the elasticity constant of the spring, k, and the mass of the body attached to the oscillator, m. This study, seeks to realize two objectives: 

1. To compare the estimation the epistemic parameter k between the samplers;

2. To compare the model updating results obtained through the use of MCMC, TMCMC, and SMC.

### 2-Dimensional Eigenvalue Problem:

This tutorial presents a 2-by-2 square [matrix](http://latex.codecogs.com/svg.latex?%5Cbegin%7Bpmatrix%7D%0D%0A%7B%5Ctheta_1%7D%2B%7B%5Ctheta_2%7D%26-%7B%5Ctheta_2%7D%5C%5C-%7B%5Ctheta_2%7D%26%7B%5Ctheta_2%7D%5C%5C%0D%0A%5Cend%7Bpmatrix%7D) in which there exists two distinct real eignvalue solutions. The matrix elements here are defined by two epistemic parameters: Theta 1 and Theta 2. This tutorial seeks to achieve three objectives:

1. To observe the performance of each of the advanced Monte Carlo samplers in obtaining samples from a 2-dimensional, bi-modal posterior distribution;

2. To estimate the solutions to the epistemic parameters: Theta 1 and Theta 2;

3. To compare the model updating results obtained through the use of MCMC, TMCMC, and SMC.

### 18-Dimensional DLR-AIRMOD Problem:

This tutorial is based on the DLR-AIRMOD problem by [Govers et. al (2014)](https://doi.org/10.1016/j.ymssp.2014.06.003). In this problem, a DLR-AIRMOD test structure is being studied and its physical structure is being defined by 18 physical parameters. Measurements are obtained from this structure in the form of its response eigen-frequencies. In total, 14 of such active eigen-frequencies have been identified. More information to the 18 parameters and 14 active eigen-frequencies can be found in the aforementioned literature. For this problem, the 14 active eigen-frequencies are mathematically defined, for the given set of 18 input parameters, by of an Artificial Neural Network surrogate model proposed by [Patelli et. al (2017)](https://doi.org/10.1016/j.proeng.2017.09.221). This surrogate model is to be updated through Bayesian Model Updating. This tutorial seeks to achieve four objectives:

1. To observe the performance of each of the advanced Monte Carlo samplers in obtaining samples from a high-dimensional posterior distribution;

2. To asses the robustness of each of the sampling algorithm;

3. To estimate the 18 epistemic input parameters;

4. To compare the model updating results obtained through the use of MCMC, TMCMC, and SMC.

Note: To run this tutorial on MATLAB, please ensure that you have installed OpenCOSSAN as well as the Fast Articifial Neural Network (FANN) Library. Instructions to install and use FANN Library can be found [HERE](https://cossan.co.uk/wiki/index.php/Installation_of_FANN_library).


## Citation

At the time of writing, the journal article to this tutorial problems presented in this repository is still under review. A pre-print version of the journal paper is made available in this repository.
