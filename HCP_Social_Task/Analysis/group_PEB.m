% Use MATLAB exports to fill in DCM objects for SPM PEB

% Subjects and conditions

load("HCP_Social_Task/Analysis/MCMC_diagnostics.mat");
conv_table = struct2table(diagnostics);

% --- Loop through subjects and conditions ---
for s = 1:height(conv_table)

    subID = conv_table.subject(s);
    cond  = conv_table.condition{s};
        
        % --- Load the R-generated posterior ---
        Rfile = sprintf('HCP_Social_Task/Output/sub_%d_%s.mat', subID, cond);
        load(Rfile);  
        
        % --- Build DCM structure ---
        DCM = struct();

        Ep = struct();
        Ep.A = A;
        Ep.B = cat(3, B1, B2);  % 3D B matrices
        Ep.C = C;

        DCM.Ep = Ep;

        % --- Set prior expectations (zeros) ---
        DCM.M.pE.A = zeros(size(A));
        DCM.M.pE.B = zeros(size(Ep.B));
        DCM.M.pE.C = zeros(size(C));

        % --- Set prior covariance mask ---
        pC_mask = zeros(size(Cp));   

        % Set the tighter variances for diagonal param as in MCMC
        idx = [1 4 9 12];
        pC_mask(idx + (idx-1)*size(pC_mask,1)) = 0.015625;

        % Set the prior var for off-diag param and C to be 1 as in MCMC
        idx = [2 3 10 11 13 14];
        pC_mask(idx + (idx-1)*size(pC_mask,1)) = 1;
        
        % Assign to DCM
        DCM.M.pC = pC_mask;

        % --- Set posterior covariance ---
        Cp = spdiags(diag(Cp), 0, size(Cp,1), size(Cp,2));
        DCM.Cp = Cp;          % sparse matrix

        % --- Free energy (ELBO) ---
        DCM.F = DCM_F;

        % --- Save DCM ---
        outFile = sprintf('HCP_Social_Task/Output/DCM_sub_%d_%s.mat', subID, cond);
        save(outFile, 'DCM', 'Ep', 'Cp', '-v7.3');
        
        fprintf('Saved %s\n', outFile);
end


%% =================
% Run PEB for each condition and export
clear; clc;
addpath /storage/work/krf5429/spm   

conditions = {'phaseLR_maskL', 'phaseLR_maskR', 'phaseRL_maskL', 'phaseRL_maskR'};

% Loop through each condition
for i = 1:numel(conditions)
    
    cond = conditions{i};
    fprintf('\n=== Running PEB for %s ===\n', cond);
    
    % 1. Select matching DCMs
    path = fullfile(pwd, 'HCP_Social_Task/Output');
    pattern = sprintf('^DCM_.*_%s\\.mat$', cond);
    DCM_files = spm_select('FPList', path, pattern);
    
    if isempty(DCM_files)
        warning('No DCMs found for pattern: %s', pattern);
        continue;
    end
    
    % 2. Design matrix
    M = struct();
    M.X = ones(size(DCM_files, 1), 1);

    % Load CSV
    data = readtable('HCP_Social_Task/HCP_YA_subjects.csv'); 
    
    % Optional: make sure 'Subject' is numeric
    if iscell(data.Subject)
        data.Subject = str2double(data.Subject);
    end
    
    N = height(data);  % number of subjects
    
    % Encode gender (assuming 1=male, 2=female, adjust if needed)
    gender = zeros(N,1);  % pre-allocate
    for k = 1:N
        if strcmp(data.Gender{k}, 'F')   % female = 1
            gender(k) = 1;
        elseif strcmp(data.Gender{k}, 'M')  % male = 0
            gender(k) = 0;
        else
            error('Unknown gender: %s', data.Gender{k});
        end
    end
    
    % Encode age ranges as ordinal
    % Define mapping
    age_map = containers.Map({'22-25','26-30','31-35','36+'}, 1:4);
    age_code = zeros(N,1);
    
    for j = 1:N
        age_code(j) = age_map(data.Age{j});
    end
    
    % Optional: z-score age
    age_z = zscore(age_code);
    
    %% Z-score PMAT
    PMAT_z = zscore(data.PMAT24_A_CR);
  
    % Include intercept + covariates
    M.X = [ones(N,1), gender, age_z, PMAT_z];
    
    % Column names
    M.Xnames = {'Intercept','Gender','Age','PMAT'};
    
    % 3. Parameter fields to include
    fields = {'A', 'B', 'C'};
    
    % 4. Estimate PEB
    PEB = spm_dcm_peb(DCM_files, M, fields);
    
    % 5. Save the full PEB result
    outname = sprintf('HCP_Social_Task/Output/PEB_%s.mat', cond);
    save(outname, 'PEB');
    
    fprintf('Saved: %s\n', outname);
end

fprintf('\nAll PEBs complete!\n');

%% ===
clear; clc;

% List of PEB files
PEB_files = {'HCP_Social_Task/Output/PEB_phaseLR_maskL.mat', 'HCP_Social_Task/Output/PEB_phaseLR_maskR.mat', ...
             'HCP_Social_Task/Output/PEB_phaseRL_maskL.mat', 'HCP_Social_Task/Output/PEB_phaseRL_maskR.mat'};

conditions = {'phaseLR_maskL', 'phaseLR_maskR', 'phaseRL_maskL', 'phaseRL_maskR'};

% Initialize master table
allPEB = table();

for f = 1:length(PEB_files)
    
    % Load PEB
    load(PEB_files{f}, 'PEB');
    
    % Number of parameters and covariates
    [nParams, nCov] = size(PEB.Ep);
    
    % Names of covariates (intercept + others)
    covNames = {'Intercept', 'Gender', 'Age', 'PMAT'};
    if nCov ~= length(covNames)
        covNames = strcat('Cov', string(1:nCov)); % fallback
    end
    
    % Extract standard deviations for each parameter & covariate
    % Cp is the posterior covariance of the vectorized Ep
    SD = sqrt(diag(PEB.Cp));
    
    % Cp corresponds to vec(Ep), so reshape into nParams x nCov
    SDmat = reshape(SD, nParams, nCov);
    
    % Make table for this condition
    T = table(PEB.Pnames(:), ...
              PEB.Ep(:,1), PEB.Ep(:,2), PEB.Ep(:,3), PEB.Ep(:,4), ...
              SDmat(:,1), SDmat(:,2), SDmat(:,3), SDmat(:,4), ...
              'VariableNames', ['Parameter', ...
                                strcat('Ep_', covNames), ...
                                strcat('SD_', covNames)]);
    
    % Add a column for condition
    T.Condition = repmat(conditions(f), nParams, 1);
    
    % Append to master table
    allPEB = [allPEB; T];
    
end

% Save table
writetable(allPEB, 'HCP_Social_Task/Analysis/Results/MCMC_PEB_Summary.csv');





