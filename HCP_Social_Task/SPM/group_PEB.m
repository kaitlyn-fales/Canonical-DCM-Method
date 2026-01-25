%==========================================================================
% Run 4 separate group PEBs (A, B, C parameters)
% No model reduction or averaging
%==========================================================================

clear; clc;
addpath /storage/work/krf5429/spm 


conditions = {'phaseLR_maskL', 'phaseLR_maskR', 'phaseRL_maskL', 'phaseRL_maskR'};

% Loop through each condition
for i = 1:numel(conditions)
    
    cond = conditions{i};
    fprintf('\n=== Running PEB for %s ===\n', cond);
    
    % 1. Select matching DCMs
    path = fullfile(pwd, 'HCP_Social_Task/SPM/Output');
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
    outname = sprintf('HCP_Social_Task/SPM/Output/PEB_%s.mat', cond);
    save(outname, 'PEB');
    
    fprintf('Saved: %s\n', outname);
end

fprintf('\nAll PEBs complete!\n');

%% ===
% List of PEB files
PEB_files = {'HCP_Social_Task/SPM/Output/PEB_phaseLR_maskL.mat', 'HCP_Social_Task/SPM/Output/PEB_phaseLR_maskR.mat', ...
             'HCP_Social_Task/SPM/Output/PEB_phaseRL_maskL.mat', 'HCP_Social_Task/SPM/Output/PEB_phaseRL_maskR.mat'};

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
writetable(allPEB, 'HCP_Social_Task/SPM/SPM_PEB_Summary.csv');



