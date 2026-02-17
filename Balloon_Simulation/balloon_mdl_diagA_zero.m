%% =============================
addpath /storage/work/krf5429/spm       

%% =============================
% General scanning specifications
n_scans = 150;
n_regions = 2;
n_inputs = 2;

% Hypothesis connectivity indices
% A: target, source
A_idxs = [2,1;
          1,2;
          1,1;
          2,2];

% B: input, target, source
B_idxs = [2,1,2];

% C: region, input
C_idxs = [1,1];

% Initialize DCM structure
DCM = struct();

% Start to populate with things that don't change between subjects and
% sessions
DCM.Y.dt = 2;          % TR (s)
DCM.Y.X0 = [];         % no confounds
DCM.v    = n_scans;    % number of scans
DCM.n    = n_regions;  % number of regions

% Initialize connectivity matrices
DCM.a = zeros(n_regions);                    % intrinsic
DCM.b = zeros(n_regions,n_regions,n_inputs); % modulatory
DCM.c = zeros(n_regions,n_inputs);           % driving

% Populate a-matrix (target, source)
for i = 1:size(A_idxs,1)
    target = A_idxs(i,1);
    source = A_idxs(i,2);
    DCM.a(target, source) = 1;
end

% Populate b-matrix (input, target, source)
for i = 1:size(B_idxs,1)
    input  = B_idxs(i,1);
    target = B_idxs(i,2);
    source = B_idxs(i,3);
    DCM.b(target, source, input) = 1;
end

% Populate c-matrix (region, input)
for i = 1:size(C_idxs,1)
    region = C_idxs(i,1);
    input  = C_idxs(i,2);
    DCM.c(region, input) = 1;
end

% Add other required elements
DCM.delays(1:n_regions) = 0; % no delays
DCM.TE     = 0.04;           % echo time

% Task-based BOLD options
DCM.options.nonlinear  = 0;       % bilinear
DCM.options.two_state  = 0;       % single-state
DCM.options.stochastic = 0;       % deterministic
DCM.options.centre     = 0;       % use U exactly as given
DCM.options.induced    = 0;       % no induced responses
DCM.options.analysis   = 'time';  % time-domain BOLD DCM
DCM.options.nograph    = 0;       % display output plots

% Define region and input names
region_names = {'V1', 'V2'};              % regions
input_names  = {'U1', 'U2'};   % experimental inputs

%% =============================
% --- Build the filename dynamically ---
filename = 'Balloon_Simulation/Data/balloon_sim_data_diagA_zero.mat';

% --- Load the file ---
load(filename);

% --- Time series data
DCM.Y.y  = Y;  

% --- xY structure
for r = 1:n_regions
    DCM.xY(r).name  = region_names{r};     % region name
    DCM.xY(r).u     = Y(:,r);              % time series
    DCM.xY(r).DT    = DCM.Y.dt;            % TR
    DCM.xY(r).XYZ   = [0 0 0];             % dummy coordinates
    DCM.xY(r).def   = 'manual';            % required
    DCM.xY(r).Sname = 'Session_1';         % dummy session
end

% --- Experimental inputs
DCM.U.u    = U;
DCM.U.name = input_names;                  % input names
DCM.U.dt   = DCM.Y.dt;

% Save and estimate
output   = 'Balloon_Simulation/Output/balloon_diagA_zero.mat';

% --- Save the DCM structure ---
save(output, 'DCM');

% --- Estimate the DCM ---
spm_dcm_estimate(output);

%% ===========
clear; clc;

% Load fitted DCM back in
load('Balloon_Simulation/Output/balloon_diagA_zero.mat');

% Save necessary variables for R
Cp = full(Cp);   

Ep_A = Ep.A;
Ep_B = Ep.B;
Ep_C = Ep.C;
DCM_y = DCM.y;

save('Balloon_Simulation/Output/balloon_diagA_zero_out.mat', 'Ep_A', 'Ep_B', 'Ep_C', 'Cp', 'DCM_y');


