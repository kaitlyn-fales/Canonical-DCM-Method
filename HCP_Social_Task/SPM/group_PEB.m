%==========================================================================
% Run 4 separate group PEBs (A, B, C parameters)
% No model reduction or averaging
%==========================================================================

clear; clc;
addpath /storage/work/krf5429/spm 
cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM/Output');  

% Define your condition suffixes
conditions = {'phaseLR_maskL', 'phaseLR_maskR', 'phaseRL_maskL', 'phaseRL_maskR'};

% Loop through each condition
for i = 1:numel(conditions)
    
    cond = conditions{i};
    fprintf('\n=== Running PEB for %s ===\n', cond);
    
    % 1. Select matching DCMs
    pattern = sprintf('^DCM_.*_%s\\.mat$', cond);
    DCM_files = spm_select('FPList', pwd, pattern);
    
    if isempty(DCM_files)
        warning('No DCMs found for pattern: %s', pattern);
        continue;
    end
    
    % 2. Group mean design matrix
    M = struct();
    M.X = ones(size(DCM_files, 1), 1);
    
    % 3. Parameter fields to include
    fields = {'A', 'B', 'C', 'transit', 'decay', 'epsilon'};
    
    % 4. Estimate PEB
    PEB = spm_dcm_peb(DCM_files, M, fields);
    
    % 5. Save the full PEB result
    outname = sprintf('PEB_%s.mat', cond);
    save(outname, 'PEB');
    
    fprintf('Saved: %s\n', outname);
end

fprintf('\nAll PEBs complete!\n');

%% ==========================
clear; clc;

phase = 'LR';    
mask  = 'L'; 

cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM/Output'); 

% Load one first-level DCM used in the PEB
load('DCM_sub_100307_phaseLR_maskL.mat')

% Load your group PEB
load('PEB_phaseLR_maskL.mat')

% Generate correct parameter names
cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM');
param_names = get_param_labels(DCM, PEB);

s = sqrt(diag(PEB.Cp)); 
s = s(1:10); 

% Combine into a labeled table
results = table(param_names(:), PEB.Ep(1:10,:), s, ...
                'VariableNames', {'Parameter', 'Mean', 'SD'});

% Define the R parameter order (names must match the R convention)
param_order_R = {
    'A(2,1)','A(1,2)','A(1,1)','A(2,2)', ...
    'B(2,2,1)','B(2,1,2)','B(2,1,1)','B(2,2,2)', ...
    'C(1,1)','C(2,1)'
};

% 1. Convert B parameter names in MATLAB table to R convention (input first)
param_names_mod = results.Parameter;  % copy MATLAB names
for r = 1:height(results)
    if startsWith(param_names_mod{r}, 'B(')
        nums = sscanf(param_names_mod{r}, 'B(%d,%d,%d)');  % target, source, input
        target = nums(1);
        source = nums(2);
        input  = nums(3);

        % reorder for R: input, target, source
        param_names_mod{r} = sprintf('B(%d,%d,%d)', input, target, source);

        % also update table columns for clarity
        results.input(r) = input;
        results.row(r)   = target;
        results.col(r)   = source;
    end
end
results.Parameter = param_names_mod;  % overwrite names with R-style

% 2. Reorder rows to match R table
[~, idx] = ismember(param_order_R, results.Parameter);
if any(idx==0)
    warning('Some R parameters were not found in MATLAB results:');
    disp(param_order_R(idx==0))
end
results_reordered = results(idx(idx>0), :);

results_final = results_reordered(:, {'Parameter', 'Mean', 'SD'});

results_final.Mean = round(results_final.Mean, 3);
results_final.SD   = round(results_final.SD, 3);

% Display the final table
msg = sprintf('PEB for Phase %s, Mask %s', phase, mask);
disp(msg)
disp(results_final)

writetable(results_final, sprintf('SPM_PEB_phase%s_mask%s.csv', phase, mask));

hemodynamics = table(PEB.Pnames(11:14), PEB.Ep(11:14,:),...
'VariableNames',{'Parameter','Mean'});

writetable(hemodynamics, sprintf('SPM_PEB_hemodynamics_phase%s_mask%s.csv', phase, mask))

%% ======================
spm_dcm_peb_review(PEB)

%% ==========================
clear; clc;

phase = 'LR';    
mask  = 'R'; 

cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM/Output'); 

% Load one first-level DCM used in the PEB
load('DCM_sub_100307_phaseLR_maskR.mat')

% Load your group PEB
load('PEB_phaseLR_maskR.mat')

% Generate correct parameter names
cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM');
param_names = get_param_labels(DCM, PEB);

s = sqrt(diag(PEB.Cp)); 
s = s(1:10); 

% Combine into a labeled table
results = table(param_names(:), PEB.Ep(1:10,:), s, ...
                'VariableNames', {'Parameter', 'Mean', 'SD'});

% Define the R parameter order (names must match the R convention)
param_order_R = {
    'A(2,1)','A(1,2)','A(1,1)','A(2,2)', ...
    'B(2,2,1)','B(2,1,2)','B(2,1,1)','B(2,2,2)', ...
    'C(1,1)','C(2,1)'
};

% 1. Convert B parameter names in MATLAB table to R convention (input first)
param_names_mod = results.Parameter;  % copy MATLAB names
for r = 1:height(results)
    if startsWith(param_names_mod{r}, 'B(')
        nums = sscanf(param_names_mod{r}, 'B(%d,%d,%d)');  % target, source, input
        target = nums(1);
        source = nums(2);
        input  = nums(3);

        % reorder for R: input, target, source
        param_names_mod{r} = sprintf('B(%d,%d,%d)', input, target, source);

        % also update table columns for clarity
        results.input(r) = input;
        results.row(r)   = target;
        results.col(r)   = source;
    end
end
results.Parameter = param_names_mod;  % overwrite names with R-style

% 2. Reorder rows to match R table
[~, idx] = ismember(param_order_R, results.Parameter);
if any(idx==0)
    warning('Some R parameters were not found in MATLAB results:');
    disp(param_order_R(idx==0))
end
results_reordered = results(idx(idx>0), :);

results_final = results_reordered(:, {'Parameter', 'Mean', 'SD'});

results_final.Mean = round(results_final.Mean, 3);
results_final.SD   = round(results_final.SD, 3);

% Display the final table
msg = sprintf('PEB for Phase %s, Mask %s', phase, mask);
disp(msg)
disp(results_final)

writetable(results_final, sprintf('SPM_PEB_phase%s_mask%s.csv', phase, mask));

hemodynamics = table(PEB.Pnames(11:14), PEB.Ep(11:14,:),...
'VariableNames',{'Parameter','Mean'});

writetable(hemodynamics, sprintf('SPM_PEB_hemodynamics_phase%s_mask%s.csv', phase, mask))

%% ======================
spm_dcm_peb_review(PEB)

%% ==========================
clear; clc;

phase = 'RL';    
mask  = 'L'; 

cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM/Output'); 

% Load one first-level DCM used in the PEB
load('DCM_sub_100307_phaseRL_maskL.mat')

% Load your group PEB
load('PEB_phaseRL_maskL.mat')

% Generate correct parameter names
cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM');
param_names = get_param_labels(DCM, PEB);

s = sqrt(diag(PEB.Cp)); 
s = s(1:10); 

% Combine into a labeled table
results = table(param_names(:), PEB.Ep(1:10,:), s, ...
                'VariableNames', {'Parameter', 'Mean', 'SD'});

% Define the R parameter order (names must match the R convention)
param_order_R = {
    'A(2,1)','A(1,2)','A(1,1)','A(2,2)', ...
    'B(2,2,1)','B(2,1,2)','B(2,1,1)','B(2,2,2)', ...
    'C(1,1)','C(2,1)'
};

% 1. Convert B parameter names in MATLAB table to R convention (input first)
param_names_mod = results.Parameter;  % copy MATLAB names
for r = 1:height(results)
    if startsWith(param_names_mod{r}, 'B(')
        nums = sscanf(param_names_mod{r}, 'B(%d,%d,%d)');  % target, source, input
        target = nums(1);
        source = nums(2);
        input  = nums(3);

        % reorder for R: input, target, source
        param_names_mod{r} = sprintf('B(%d,%d,%d)', input, target, source);

        % also update table columns for clarity
        results.input(r) = input;
        results.row(r)   = target;
        results.col(r)   = source;
    end
end
results.Parameter = param_names_mod;  % overwrite names with R-style

% 2. Reorder rows to match R table
[~, idx] = ismember(param_order_R, results.Parameter);
if any(idx==0)
    warning('Some R parameters were not found in MATLAB results:');
    disp(param_order_R(idx==0))
end
results_reordered = results(idx(idx>0), :);

results_final = results_reordered(:, {'Parameter', 'Mean', 'SD'});

results_final.Mean = round(results_final.Mean, 3);
results_final.SD   = round(results_final.SD, 3);

% Display the final table
msg = sprintf('PEB for Phase %s, Mask %s', phase, mask);
disp(msg)
disp(results_final)

writetable(results_final, sprintf('SPM_PEB_phase%s_mask%s.csv', phase, mask));

hemodynamics = table(PEB.Pnames(11:14), PEB.Ep(11:14,:),...
'VariableNames',{'Parameter','Mean'});

writetable(hemodynamics, sprintf('SPM_PEB_hemodynamics_phase%s_mask%s.csv', phase, mask))

%% ======================
spm_dcm_peb_review(PEB)

%% ==========================
clear; clc;

phase = 'RL';    
mask  = 'R'; 

cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM/Output'); 

% Load one first-level DCM used in the PEB
load('DCM_sub_100307_phaseRL_maskR.mat')

% Load your group PEB
load('PEB_phaseRL_maskR.mat')

% Generate correct parameter names
cd('/storage/work/krf5429/Canonical-DCM-Method/HCP_Social_Task/SPM');
param_names = get_param_labels(DCM, PEB);

s = sqrt(diag(PEB.Cp)); 
s = s(1:10); 

% Combine into a labeled table
results = table(param_names(:), PEB.Ep(1:10,:), s, ...
                'VariableNames', {'Parameter', 'Mean', 'SD'});

% Define the R parameter order (names must match the R convention)
param_order_R = {
    'A(2,1)','A(1,2)','A(1,1)','A(2,2)', ...
    'B(2,2,1)','B(2,1,2)','B(2,1,1)','B(2,2,2)', ...
    'C(1,1)','C(2,1)'
};

% 1. Convert B parameter names in MATLAB table to R convention (input first)
param_names_mod = results.Parameter;  % copy MATLAB names
for r = 1:height(results)
    if startsWith(param_names_mod{r}, 'B(')
        nums = sscanf(param_names_mod{r}, 'B(%d,%d,%d)');  % target, source, input
        target = nums(1);
        source = nums(2);
        input  = nums(3);

        % reorder for R: input, target, source
        param_names_mod{r} = sprintf('B(%d,%d,%d)', input, target, source);

        % also update table columns for clarity
        results.input(r) = input;
        results.row(r)   = target;
        results.col(r)   = source;
    end
end
results.Parameter = param_names_mod;  % overwrite names with R-style

% 2. Reorder rows to match R table
[~, idx] = ismember(param_order_R, results.Parameter);
if any(idx==0)
    warning('Some R parameters were not found in MATLAB results:');
    disp(param_order_R(idx==0))
end
results_reordered = results(idx(idx>0), :);

results_final = results_reordered(:, {'Parameter', 'Mean', 'SD'});

results_final.Mean = round(results_final.Mean, 3);
results_final.SD   = round(results_final.SD, 3);

% Display the final table
msg = sprintf('PEB for Phase %s, Mask %s', phase, mask);
disp(msg)
disp(results_final)

writetable(results_final, sprintf('SPM_PEB_phase%s_mask%s.csv', phase, mask));

hemodynamics = table(PEB.Pnames(11:14), PEB.Ep(11:14,:),...
'VariableNames',{'Parameter','Mean'});

writetable(hemodynamics, sprintf('SPM_PEB_hemodynamics_phase%s_mask%s.csv', phase, mask))

%% ======================
spm_dcm_peb_review(PEB)
