

close all
clear variables

%-----------------------------------------------------------
%-------------- DCM model specification --------------------

load("y_obs_SPM.mat");
load("u_full_SPM.mat");

%--- Basic settings
f_fname = @f_DCMwHRF;
g_fname = @g_HRF3;
TR = 3.22;                    % sampling period (in sec)
n_t = 360;                    % number of time samples             
microDT = dt;                % micro-time resolution (in sec)
homogeneous = 0;              % params of g(x) homogeneous accross regions
reduced_f = 0;                % fix some HRF params
lin = 1;                      % linearized variant of HRF Balloon model
stochastic = 0;               % flag for stochastic DCM inversion
alpha = Inf;%1e2/TR;               % state noise precision
sigma = 1e0;                  % measurement noise precision

y = transpose(y_obs);
u = transpose(u_full);
nu = size(u,1);

%--- DCM structure
% invariant effective connectivity
A = [0 1 0
     1 0 1
     0 1 0];
nreg = size(A,1);
% modulatory effects
B{1} = zeros(nreg,nreg);
B{2} = [0 0 0
        1 0 0
        0 0 0];
B{3} = [0 0 0
        1 0 0
        0 0 0];
% input-state coupling
C = [1 0 0
     0 0 0
     0 0 0];
% gating (nonlinear) effects
D{1} = zeros(nreg,nreg);
D{2} = zeros(nreg,nreg);
D{3} = zeros(nreg,nreg);

%--- Build options and dim structures
options = prepare_fullDCM(A,B,C,D,TR,microDT,homogeneous);
options.priors = getPriors(nreg,n_t,options,reduced_f,stochastic);
options.microU = 1;
options.backwardLag = ceil(16/TR);  % 16 secs effective backward lag
options.inF.linearized = lin;
dim.n_theta = options.inF.ind5(end);
if options.inG.homogeneous
    dim.n_phi = 2;
else
    dim.n_phi = 2*nreg;
end  
dim.n = 5*nreg;


%--- Call inversion routine
[p2,o2] = VBA_NLStateSpaceModel(y,u,f_fname,g_fname,dim,options);
VBA_spm_dcm_explore(vba2dcm(p2,o2,[],TR));

% Save
save VBA_results.mat p2
