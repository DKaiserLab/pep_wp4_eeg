function ds = makeChunkForIncompleteDS(ds, cfg, iSub)

% this checks if the data set is complete and if not takes care of this
% Note that for every participant with missing data an individual solution
% must be found

% kick out first block (for subject 114) - connection problem in
% block 1
if cfg.subNums(iSub) == 114
    % get number of missing trials
    missing_trials = cfg.nTrials*cfg.nBlocks - height(ds.samples);

    % add fake trials
    ds.samples(1+missing_trials:end+missing_trials, :) = ds.samples;
    ds.sa.trialinfo(1+missing_trials:end+missing_trials) = ds.sa.trialinfo;
    ds.sa.trialinfo(1:10) = nan;
    ds.sa.trialinfo(1:10) = setdiff(unique(1:cfg.nTrials/2), unique(ds.sa.trialinfo(1:cfg.nTrials/2)'));

    % make chunks
    nch=cfg.nBlocks/cfg.nTrials2average;
    ds.sa.chunks=(1:length(ds.sa.trialinfo))';
    ds.sa.targets=ds.sa.trialinfo;
    ds.sa.chunks=cosmo_chunkize(ds, nch);

    if cfg.nTrials2average == 1

        % remove first block
        includeTrials = ones(1, cfg.nTrials*cfg.nBlocks);
        includeTrials(1:cfg.nTrials/2) = 0;
        ds = cosmo_slice(ds, logical(includeTrials'));

    elseif cfg.nTrials2average == 2


        % replace fake trials by copy of trial they get averaged with
        ds.samples(1:cfg.nTrials/2, :) = ds.samples(cfg.nTrials*nch+1:cfg.nTrials*nch+cfg.nTrials/2, :);
        ds.sa.targets(1:cfg.nTrials/2) = ds.sa.targets(cfg.nTrials*nch+1:cfg.nTrials*nch+cfg.nTrials/2);
        ds.sa.trialinfo(1:cfg.nTrials/2) = ds.sa.trialinfo(cfg.nTrials*nch+1:cfg.nTrials*nch+cfg.nTrials/2);

        % averaging trials
        ds = cosmo_average_samples(ds, 'targets', ds.sa.targets, 'chunks', ds.sa.chunks);

    else
        warning('Find solution for this subject and this trial averaging condition')

    end


else
    warning('Find solution for this subject')
end
end