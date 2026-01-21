function param_names = get_param_labels(DCM, PEB)
% Returns parameter names for PEB.Ep using DCM.a, DCM.b, DCM.c
% Works for partial B/C and multi-input 3D B matrices

param_names = {};
counter = 0;
max_params = length(PEB.Ep);

% --- A
for i = 1:size(DCM.M.pE.A,1)
    for j = 1:size(DCM.M.pE.A,2)
        if DCM.a(i,j) == 1
            counter = counter + 1;
            param_names{counter} = sprintf('A(%d,%d)', i,j);
            if counter >= max_params, return; end
        end
    end
end

% --- B (3D: target x source x input)
Bmask = DCM.b;
for k = 1:size(Bmask,3)       % input
    for j = 1:size(Bmask,2)   % source
        for i = 1:size(Bmask,1) % target
            if Bmask(i,j,k) == 1
                counter = counter + 1;
                param_names{counter} = sprintf('B(%d,%d,%d)', i,j,k);
                if counter >= max_params, return; end
            end
        end
    end
end

% --- C
for i = 1:size(DCM.M.pE.C,1)
    for j = 1:size(DCM.M.pE.C,2)
        if DCM.c(i,j) == 1
            counter = counter + 1;
            param_names{counter} = sprintf('C(%d,%d)', i,j);
            if counter >= max_params, return; end
        end
    end
end

end
