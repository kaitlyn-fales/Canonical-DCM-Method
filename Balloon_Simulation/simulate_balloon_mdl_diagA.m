%% =============================
% Balloon model with diagonal parameters on A only and no delays

clear;clc;

% General scanning specifications
n_scans = 150;
n_regions = 2;
n_inputs = 2;

% Initialize DCM structure
DCM = struct();

DCM.Y.dt = 2;          % TR (s)
DCM.Y.X0 = [];         % no confounds
DCM.v    = n_scans;    % number of scans
DCM.n    = n_regions;  % number of regions

% Initialize connectivity matrices for hypothesis
DCM.a = ones(n_regions);                     % intrinsic - all connections
DCM.b = zeros(n_regions,n_regions,n_inputs); % modulatory
DCM.c = zeros(n_regions,n_inputs);           % driving
DCM.d = zeros(n_regions,n_regions,0);          % non-linear

% Add in non-zero B and C
DCM.b(1,2,2) = 1;
DCM.c(1,1) = 1;

% A matrix
A = [-0.1  0.3;
      0.4  0.15];

% B matrices (array)
B = zeros(2, 2, 2);  

B(:,:,1) = [0   0;
            0   0];

B(:,:,2) = [0  -0.2;
            0   0];

% C matrix
C = [0.7  0;
     0.0  0];

% Add parameters into DCM structure
DCM.Ep.A = A;
DCM.Ep.B = B;
DCM.Ep.C = C;
DCM.Ep.D = DCM.d;

% Hemodynamic parameters
DCM.Ep.transit = [-0.2; -0.3];   % tau
DCM.Ep.decay   = -0.2;           % kappa
DCM.Ep.epsilon = 0.15;           % epsilon

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
input_names  = {'Input1', 'Input2'};   % experimental inputs

% --- Time series data - placeholder
DCM.Y.y    = [];

% --- xY structure
for r = 1:n_regions
    DCM.xY(r).name  = region_names{r};     % region name
    DCM.xY(r).u     = [];              % time series
    DCM.xY(r).DT    = DCM.Y.dt;            % TR
    DCM.xY(r).XYZ   = [0 0 0];             % dummy coordinates
    DCM.xY(r).def   = 'manual';            % required
    DCM.xY(r).Sname = 'Session_1';         % dummy session
end

% --- Experimental inputs
load("Balloon_Simulation/Data/U_mat.mat");
DCM.U.u    = U;
DCM.U.name = input_names;                  % input names
DCM.U.dt   = DCM.Y.dt;

% Simulate data
[Y, ~, ~] = spm_dcm_generate(DCM, inf, false);

% Save output simulated data and U matrix for R to then add noise to
y_signal = Y.y;  % extract BOLD signal
save('Balloon_Simulation/Data/balloon_sim_signal_diagA_zero.mat', 'y_signal', 'U');

%% ==============================
% Balloon model with diagonal parameters on A only and 1s delays (TR/2)

clear; clc;

% General scanning specifications
n_scans = 150;
n_regions = 2;
n_inputs = 2;

% Initialize DCM structure
DCM = struct();

DCM.Y.dt = 2;          % TR (s)
DCM.Y.X0 = [];         % no confounds
DCM.v    = n_scans;    % number of scans
DCM.n    = n_regions;  % number of regions

% Initialize connectivity matrices for hypothesis
DCM.a = ones(n_regions);                     % intrinsic - all connections
DCM.b = zeros(n_regions,n_regions,n_inputs); % modulatory
DCM.c = zeros(n_regions,n_inputs);           % driving
DCM.d = zeros(n_regions,n_regions,0);          % non-linear

% Add in non-zero B and C
DCM.b(1,2,2) = 1;
DCM.c(1,1) = 1;

% A matrix
A = [-0.1  0.3;
      0.4  0.15];

% B matrices (array)
B = zeros(2, 2, 2);  

B(:,:,1) = [0   0;
            0   0];

B(:,:,2) = [0  -0.2;
            0   0];

% C matrix
C = [0.7  0;
     0.0  0];

% Add parameters into DCM structure
DCM.Ep.A = A;
DCM.Ep.B = B;
DCM.Ep.C = C;
DCM.Ep.D = DCM.d;

% Hemodynamic parameters
DCM.Ep.transit = [-0.2; -0.3];   % tau
DCM.Ep.decay   = -0.2;           % kappa
DCM.Ep.epsilon = 0.15;           % epsilon

% Add other required elements
DCM.delays(1:n_regions) = DCM.Y.dt / 2; % 1s delays
DCM.TE     = 0.04;                      % echo time

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
input_names  = {'Input1', 'Input2'};   % experimental inputs

% --- Time series data - placeholder
DCM.Y.y    = [];

% --- xY structure
for r = 1:n_regions
    DCM.xY(r).name  = region_names{r};     % region name
    DCM.xY(r).u     = [];              % time series
    DCM.xY(r).DT    = DCM.Y.dt;            % TR
    DCM.xY(r).XYZ   = [0 0 0];             % dummy coordinates
    DCM.xY(r).def   = 'manual';            % required
    DCM.xY(r).Sname = 'Session_1';         % dummy session
end

% --- Experimental inputs
load("Balloon_Simulation/Data/U_mat.mat");
DCM.U.u    = U;
DCM.U.name = input_names;                  % input names
DCM.U.dt   = DCM.Y.dt;

% Simulate data
[Y, ~, ~] = spm_dcm_generate(DCM, inf, false);

% Save output simulated data and U matrix for R to then add noise to
y_signal = Y.y;  % extract BOLD signal
save('Balloon_Simulation/Data/balloon_sim_signal_diagA_nonzero.mat', 'y_signal', 'U');
