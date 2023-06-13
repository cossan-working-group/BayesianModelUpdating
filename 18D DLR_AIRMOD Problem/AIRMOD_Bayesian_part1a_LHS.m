Spath = fileparts(which('AIRMOD_Bayesian_part1a_LHS.m'));
OpenCossan.setWorkingPath(fullfile(Spath,'workfolder'));

%% Model definition
%% Connector definition

Xinj = Injector('Sscanfilepath',Spath,...
    'Sscanfilename','PARAM.DAT.cossan',...
    'Sfile','PARAM.DAT');

% Number of extracted modes:
Nmodes = 30; 

% Calculating eigenfrequencies of each extracted mode:
Xresponse = Response(...
    'Sname','eigenfrequencies',...
    'Clookoutfor',{'         R E A L   E I G E N V A L U E S'},...
    'Nrow',6, 'Ncol',68,...
    'Nrepeat',Nmodes,...
    'Sformat','%12e');

Xext = Extractor('Sfile','AIRMOD103.f06',...
    'CXresponse',{Xresponse});

Xconnector = Connector('Stype','nastran',...
    'Ssolverbinary','/usr/software/Nastran/20122/bin/nast20122',... '/Apps/msc/MSC_Nastran/20131/bin/nast20131',...
    'Sexecmd','%Ssolverbinary %Smaininputfile %Sexeflags',...
    'Sexeflags','scr=yes news=no bat=no old=no',...
    'SerrorFileExtension','f06',...
    'SerrorString','FATAL',...
    'SmainInputPath',Spath,...
    'SmainInputFile','AIRMOD103.bdf',...
    'CSadditionalFiles',{'AIRMOD_SOLID_MODEL_ONLY.BDF'},...
    'LkeepSimulationFiles',false);

Xconnector = Xconnector.add(Xinj);
Xconnector = Xconnector.add(Xext);
%% Post-processing MIO
% This script splits the eigenfrequencies, saved in a dataseries by the
% connector, into individual scalar outputs:

CeigfNames = cell(1,Nmodes);
for i = 1:Nmodes
    CeigfNames{i} = ['freq' num2str(i)];
end

Xmio = Mio('Spath',Spath,'Sfile','mio_post.m',...
    'CinputNames',{'eigenfrequencies'},...
    'Coutputnames',CeigfNames,...
    'Lfunction',false,'LioStructure',true,'LioMatrix',false);


%% Evaluator
Xeval = Evaluator('CXmembers',{Xconnector,Xmio},'CSmembers',{'Xconnector','Xmio'},...
    'XjobManagerInterface',JobManagerInterface('Stype','GridEngine'),...
    'CSqueues',{'all.q',''},'CShostnames',{'cossan.cfd.liv.ac.uk',''},...
    'Vconcurrent',[25 0],'LremoteInjectExtract',true);
%% Prior
% Define the prior random variables:

% We set the lowerbound to 95% under the mean and upperbound to 100% over
% the mean. The magic numbers come from:

[mean_factor,v] = unifstat(0.05,2.0); % normalized bounds
CoV = sqrt(v)/mean_factor;

% Support stiffness:
theta01 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.8e3,'cov',CoV); TinputDefault.theta01 = 1.8e3;
theta02 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.5e3,'cov',CoV); TinputDefault.theta02 = 7.5e3;
theta03 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.3e2,'cov',CoV); TinputDefault.theta03 = 1.3e2;
theta04 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.0e1,'cov',CoV); TinputDefault.theta04 = 7.0e1;
theta05 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.0e1,'cov',CoV); TinputDefault.theta05 = 7.0e1;

% Joint stiffness:
theta06 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.0e7,'cov',CoV); TinputDefault.theta06 = 1.0e7;
theta07 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.0e9,'cov',CoV); TinputDefault.theta07 = 1.0e9;
theta14 = RandomVariable('Sdistribution','uniform','mean',mean_factor*2.0e7,'cov',CoV); TinputDefault.theta14 = 2.0e7;
theta15 = RandomVariable('Sdistribution','uniform','mean',mean_factor*2.0e7,'cov',CoV); TinputDefault.theta15 = 2.0e7;
theta16 = RandomVariable('Sdistribution','uniform','mean',mean_factor*7.0e6,'cov',CoV); TinputDefault.theta16 = 7.0e6;
theta17 = RandomVariable('Sdistribution','uniform','mean',mean_factor*5.0e7,'cov',CoV); TinputDefault.theta17 = 5.0e7;
theta18 = RandomVariable('Sdistribution','uniform','mean',mean_factor*5.0e7,'cov',CoV); TinputDefault.theta18 = 5.0e7;

% Masses (sensor cables, screws and glue):
theta08 = RandomVariable('Sdistribution','uniform','mean',mean_factor*2.0e-1,'cov',CoV); TinputDefault.theta08 = 2.0e-1;
theta09 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.86e-1,'cov',CoV); TinputDefault.theta09 = 1.86e-1;
theta10 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.86e-1,'cov',CoV); TinputDefault.theta10 = 1.86e-1;
theta11 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.5e-2,'cov',CoV); TinputDefault.theta11 = 1.5e-2;
theta12 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.5e-2,'cov',CoV); TinputDefault.theta12 = 1.5e-2;
theta13 = RandomVariable('Sdistribution','uniform','mean',mean_factor*1.5e-2,'cov',CoV); TinputDefault.theta13 = 1.5e-2;

Xrvset = RandomVariableSet('Cmembers',{'theta01','theta02','theta03','theta04',...
    'theta05','theta06','theta07','theta08','theta09','theta10','theta11',...
    'theta12','theta13','theta14','theta15','theta16','theta17','theta18'});

Xinput = Input('CXmembers',{Xrvset},'CSmembers',{'Xrvset'});

%% Test model

Xmodel = Model('Xinput',Xinput,'Xevaluator',Xeval);

% Save the Test model data:
save fullmodel.mat Xmodel

%% Monte-carlo method:

% Latin Hypercube sampling:
Xlhs = LatinHypercubeSampling('Nsamples',5000,'Nbatches',10);
Xout_LHS = Xlhs.apply(Xmodel);
SbatchFolder = Xout_LHS.SbatchFolder;

%% Save the samples from Latin Hypercube sampling:

save fullmodel.mat Xmodel Xlhs SbatchFolder