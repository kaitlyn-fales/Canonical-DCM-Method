addpath /storage/work/krf5429/spm 

subjects = [ ...
100307,100408,101107,101309,101915,103111,103414,103818,105014, ...
105115,106016,108828,110411,111312,111716,113619,113922,114419, ...
115320,116524,117122,118528,118730,118932,120111,122317,122620, ...
123117,123925,124422,125525,126325,127630,127933,128127,128632, ...
129028,130013,130316,131217,131722,133019,133928,135225,135932, ...
136833,138534,139637,140925,144832,146432,147737,148335,148840, ...
149337,149539,149741,151223,151526,151627,153025,154734,156637, ...
159340,160123,161731,162733,163129,176542,178950,188347,189450, ...
190031,192540,196750,198451,199655,201111,208226,211417,211720, ...
212318,214423,221319,239944,245333,280739,298051,366446,397760, ...
414229,499566,654754,672756,751348,756055,792564,856766,857263, ...
899885];

conditions = {'phaseLR_maskL', 'phaseLR_maskR', 'phaseRL_maskL', 'phaseRL_maskR'};

for s = 1:length(subjects)
    subj = subjects(s);

    for c = 1:length(conditions)
        cond = conditions{c};

        fprintf('Processing sub_%d, condition %s...\n', subj, cond);

        dcm_file = sprintf('HCP_Social_Task/SPM/Output/DCM_sub_%d_%s.mat', subj, cond);
        load(dcm_file, 'DCM', 'Ep');

        % Hemodynamic parameters
        transit = full(Ep.transit);   % tau
        decay   = full(Ep.decay);        % kappa
        epsilon = full(Ep.epsilon);        % epsilon
        
        % U matrix
        U = DCM.U.u;
        
        % Clear environment of fitted DCM
        clear Cp DCM Ep F

        %% ==================================
        % Load in group level posterior means, restructure, and transform diagonal
        % back
        load("HCP_Social_Task/Analysis/Results/MCMC_results.mat");
        
        % Choose group level phase and mask
        results = str2double(MCMC_results.(cond)(:,2));
        
        % Reshape results into A, B, C
        A = reshape(results(1:4), 2, 2);
        
        B = zeros(2, 2, 2); 
        B(:,:,2) = reshape(results(5:8), 2, 2);
        
        C = zeros(2, 2);    
        C(:,1) = results(9:10); 
        
        % Transform diagonal of A back to raw values (invert
        % reparameterization)
        a_diag = diag(A);
        a_raw = log(a_diag ./ (-0.5));
        
        % Put back into matrices
        A(logical(eye(size(A)))) = a_raw;   

        clear MCMC_results results a_diag b_diag a_raw b_raw temp

        %% ==================================
        
        % General scanning specifications
        n_scans = 274;
        n_regions = 2;
        n_inputs = 2;
        
        % Initialize DCM structure
        DCM = struct();
        
        DCM.Y.dt = 0.72;          % TR (s)
        DCM.Y.X0 = [];         % no confounds
        DCM.v    = n_scans;    % number of scans
        DCM.n    = n_regions;  % number of regions
        DCM.U.dt = DCM.Y.dt;
        
        % Initialize connectivity matrices for hypothesis
        DCM.a = ones(n_regions);                     % intrinsic - all connections
        DCM.b(:,:,1) = zeros(n_regions,n_regions);
        DCM.b(:,:,2) = ones(n_regions,n_regions); % modulatory
        DCM.c = zeros(n_regions,n_inputs);           % driving
        DCM.d = zeros(n_regions,n_regions,0);          % non-linear
        
        % Add in non-zero C
        DCM.c(1,1) = 1;
        DCM.c(2,1) = 1;
        
        % Add parameters into DCM structure
        DCM.Ep.A = A;
        DCM.Ep.B = B;
        DCM.Ep.C = C;
        DCM.Ep.D = DCM.d;
        
        DCM.Ep.transit = transit;
        DCM.Ep.decay = decay;
        DCM.Ep.epsilon = epsilon;
        
        % TE and delays (SPM standard delay is TR/2)
        DCM.TE = 0.0331;
        DCM.delays(1:n_regions) = 0.36;
        
        % Task-based BOLD options
        DCM.options.nonlinear  = 0;       % bilinear
        DCM.options.two_state  = 0;       % single-state
        DCM.options.stochastic = 0;       % deterministic
        DCM.options.centre     = 0;       % use U exactly as given
        DCM.options.induced    = 0;       % no induced responses
        DCM.options.analysis   = 'time';  % time-domain BOLD DCM
        DCM.options.nograph    = 0;       % display output plots
        
        % Define region and input names
        region_names = {'V5', 'pSTS'};              % regions
        input_names  = {'All Motion', 'Animate'};   % experimental inputs
        
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
        DCM.U.u    = U;
        DCM.U.name = input_names;                  % input names
        DCM.U.dt   = DCM.Y.dt;

        clc;
        clearvars -except DCM U subj cond s c subjects conditions
        
        %% =====
        % Simulate data
        [Y, ~, ~] = spm_dcm_generate(DCM, inf, false);

        % Save simulated signal
        y_signal = Y.y;  
        save_path = sprintf('HCP_Social_Task/SPM/Sensitivity_Sim/MCMCpar_SPMdgm/Data/sub_%d_%s.mat', subj, cond);
        save(save_path, 'y_signal', 'U');

    end
end

