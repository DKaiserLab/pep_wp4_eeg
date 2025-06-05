% define parameters
cfg.subNums=[102 103 104 106 109 110 111 112 113 115 116 117 120 122];
cfg.makeBetweenComparison = false;
cfg.nTrials=100;

% evaluate input
if ~isfield(cfg, 'plotting'); cfg.plotting = true; end
if ~isfield(cfg, 'makeBetweenComparison'); cfg.makeBetweenComparison = false; end
makeBetweenComparison = cfg.makeBetweenComparison;

% Define classifiers
classifier = @cosmo_classify_lda;

% get existing pairwise decoding results
fileName = fullfile(pwd, '..', '..', 'derivatives', 'group_level', 'RDM',...
    'results_RDM_of_pairwise_decoding.mat');
if exist(fileName, 'file')
    load(fileName)
end

% loop through subjects
for iSub = 1:length(cfg.subNums)

    subID = sprintf('sub-%0.3d', cfg.subNums(iSub));

    % check if is part of result structure already
    subID2 = strrep(subID, '-', '');
    if exist('res', 'var')
        if isfield(res, subID2)
            disp(['RDM for ', subID, ' already exists']);
            continue
        end
    end

    % progress report
    disp(['Starting pairwise decoding for subject ',  num2str(cfg.subNums(iSub))]);

    filepath = fullfile(pwd, '..', '..', 'derivatives', ['sub-', num2str(cfg.subNums(iSub))], 'eeg', ['PEP_WP4_EEG', num2str(cfg.subNums(iSub)), '_timelock', '.mat']);
    load(filepath);

    %convert to cosmo
    ds=cosmo_meeg_dataset(timelock);
    %clear timelock

    % time range for decoding
    decoding_start = 0;
    decoding_end = 0.3;
    time_points = find(ds.a.fdim.values{2, 1} >= decoding_start &...
        ds.a.fdim.values{2, 1} <= decoding_end);

    %get time info
    res.all_time=timelock.time;
    res.included_time=time_points;
    clear timelock%new


    %% Pairwise decoding
    % Initialize RDM
    rdm = zeros(cfg.nTrials, cfg.nTrials, length(time_points));
    mean_accuracy = zeros(1, length(time_points));
    disp('')

    for tp=1:length(time_points)%1:max(ds.fa.time)
        tp_idx = time_points(tp);


        disp(char(datetime))

        if isempty(gcp('nocreate'))
            parpool(8);
        end
        nTrials = cfg.nTrials;
        parfor stim1 = 1:nTrials
            for stim2 = 1:nTrials
                if ~(stim2 > stim1)
                    continue
                end

                % check whether its a between stimulus comparison
                if ~makeBetweenComparison
                    if (stim2 > nTrials/2 && stim1 <= nTrials/2)
                        rdm(stim2, stim1, tp) = NaN;
                        continue
                    end
                end

                % Subset data for the two stimuli
                ds_stim=cosmo_slice(ds, ds.sa.trialinfo == stim1 | ds.sa.trialinfo == stim2);
                ds_stim=cosmo_slice(ds_stim, ismember(ds.fa.time,tp_idx),2);

                % Rename target
                ds_stim.sa.targets = (ds_stim.sa.trialinfo == stim1) + 1;


                nch=20;

                ds_stim.sa.chunks=[1:length(ds_stim.sa.targets)]';
                ds_stim.sa.chunks=cosmo_chunkize(ds_stim,nch);

                % Define partitions
                partitions = cosmo_nchoosek_partitioner(ds_stim, 1);

                % get predictions for each fold
                [~ ,accuracy] = cosmo_crossvalidate(ds_stim, classifier, partitions);

                % Store the accuracy
                rdm(stim2, stim1, tp) = accuracy;
            end
        end

        % make rdm symetric
        rdm(:, :, tp) = squareform(squareform(rdm(:, :, tp)));

        % take mean
        mean_accuracy(tp) = mean(squareform(rdm(:, :, tp)), 'omitnan');

    end

    % Save the RDM
    res.(subID2).rdm = rdm;
    res.(subID2).mean_accuracy = mean_accuracy;
end

% save results
outputFolder = fullfile(pwd, '..', '..', 'derivatives', 'group_level', 'RDM');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
save(fileName, 'res')

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end