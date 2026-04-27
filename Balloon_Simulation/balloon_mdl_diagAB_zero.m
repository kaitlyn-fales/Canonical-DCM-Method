clear; clc;

addpath /storage/work/krf5429/spm

n_reps = 50;

%% General scanning specifications
n_scans   = 150;
n_regions = 2;
n_inputs  = 2;

A_idxs = [2,1;
          1,2;
          1,1;
          2,2];

B_idxs = [2,1,2;
          2,2,2];

C_idxs = [1,1];

region_names = {'V1', 'V2'};
input_names  = {'U1', 'U2'};

for rep_id = 1:n_reps

    fprintf('Running replicate %03d\n', rep_id);

    %% Initialize DCM
    DCM = struct();

    DCM.Y.dt = 2;
    DCM.Y.X0 = [];
    DCM.v    = n_scans;
    DCM.n    = n_regions;

    DCM.a = zeros(n_regions);
    DCM.b = zeros(n_regions,n_regions,n_inputs);
    DCM.c = zeros(n_regions,n_inputs);

    for i = 1:size(A_idxs,1)
        target = A_idxs(i,1);
        source = A_idxs(i,2);
        DCM.a(target, source) = 1;
    end

    for i = 1:size(B_idxs,1)
        input  = B_idxs(i,1);
        target = B_idxs(i,2);
        source = B_idxs(i,3);
        DCM.b(target, source, input) = 1;
    end

    for i = 1:size(C_idxs,1)
        region = C_idxs(i,1);
        input  = C_idxs(i,2);
        DCM.c(region, input) = 1;
    end

    DCM.delays(1:n_regions) = 0;
    DCM.TE = 0.04;

    DCM.options.nonlinear  = 0;
    DCM.options.two_state  = 0;
    DCM.options.stochastic = 0;
    DCM.options.centre     = 0;
    DCM.options.induced    = 0;
    DCM.options.analysis   = 'time';
    DCM.options.nograph    = 1;

    %% Load replicate data
    filename = sprintf( ...
        'Balloon_Simulation/Data/balloon_sim_diagAB_zero_data_rep%03d.mat', ...
        rep_id ...
    );

    load(filename, 'U', 'Y');

    DCM.Y.y = Y;

    for r = 1:n_regions
        DCM.xY(r).name  = region_names{r};
        DCM.xY(r).u     = Y(:,r);
        DCM.xY(r).DT    = DCM.Y.dt;
        DCM.xY(r).XYZ   = [0 0 0];
        DCM.xY(r).def   = 'manual';
        DCM.xY(r).Sname = 'Session_1';
    end

    DCM.U.u    = U;
    DCM.U.name = input_names;
    DCM.U.dt   = DCM.Y.dt;

    %% Save and estimate
    output = sprintf( ...
        'Balloon_Simulation/Output/balloon_diagAB_zero_rep%03d.mat', ...
        rep_id ...
    );

    save(output, 'DCM');

    spm_dcm_estimate(output);

    %% Load fitted DCM and save compact output
    fitted = load(output);

    Cp = full(fitted.Cp);

    Ep_A = fitted.Ep.A;
    Ep_B = fitted.Ep.B;
    Ep_C = fitted.Ep.C;
    DCM_y = fitted.DCM.y;

    output_compact = sprintf( ...
        'Balloon_Simulation/Output/balloon_diagAB_zero_rep%03d_out.mat', ...
        rep_id ...
    );

    save(output_compact, 'Ep_A', 'Ep_B', 'Ep_C', 'Cp', 'DCM_y');

end