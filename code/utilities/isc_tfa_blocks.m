function d = isc_tfa_blocks(d, cfg)

if ~isfield(cfg, 'threshold'); cfg.threshold = 0.1; end
if ~isfield(cfg, 'nBlks'); cfg.nBlks = 20; end
if ~isfield(cfg, 'downsampleN'); cfg.downsampleN = 1; end % downsampleN = 1 -> no downsampling
if ~isfield(cfg, 'frequencies'); cfg.frequencies = {'delta', 'theta', 'alpha', 'beta', 'gamma'}; end

%% get cleaned subjects data

allTfFile = fullfile(cfg.outputPath, 'group_level', 'tfa', 'averageTfAll.mat');
if exist(allTfFile, 'file')
    load(allTfFile)
else
    allTfOut = clean_tfa_blocks(cfg);
end


for frq = 1:length(cfg.frequencies)
    frqBand = cfg.frequencies{frq};
    frqIdx = find(strcmp(allTfOut{1}.tf_all{1}.bandlabel, frqBand));

    kitMat = nan(cfg.n, cfg.n, cfg.nBlks/2);
    batMat = kitMat;
    kitMedians = nan(1, cfg.nBlks/2);
    batMedians = kitMedians;

    % loop through blocks
    for b = 1:cfg.nBlks

        if mod(b,2) == 0
            category = 'kitchen';
        else
            category = 'bathroom';
        end


        tps = downsample(allTfOut{1}.tf_all{b}.time, cfg.downsampleN);
        timecourses = nan(length(tps), cfg.n);

        % loop through subjects
        for iSub = 1:cfg.n
            subID = sprintf('sub-%0.3d', cfg.subNums(iSub));
            subID2 = strrep(subID, '-', '');

            % make a matrix with vectorized RDMs
            timecourse = squeeze(allTfOut{iSub}.tf_all{b}.powspctrm(:,:,frqIdx,:));

            % downsample
            timecourses(:, iSub) = downsample(timecourse, cfg.downsampleN);
        end

        % remove nan rows
        timecourses = timecourses(min(~isnan(timecourses), [], 2), :);

        % make IS-RDM
        [~, mat_out, ~] = make_RDM(timecourses, cfg);
        if cfg.dissimilarity
            median_mat_out = 1 - mat_out;
        else
            median_mat_out = mat_out;
        end

        % store in structure
        d.ISC.([category,'_RDM']).(frqBand).blockTfa(ceil(b/2)).RDM = mat_out;
        d.ISC.([category,'_RDM']).(frqBand).blockTfa(ceil(b/2)).color = [0, 0, 0];
        d.ISC.([category,'_RDM']).(frqBand).blockTfa(ceil(b/2)).name = (num2str(b));
        d.ISC.medianISC.(category).(frqBand).blockTfa(ceil(b/2)).medianISC = ...
            median(median_mat_out, 'all', 'omitnan');

        if mod(b,2) == 0
            kitMat(:, :, ceil(b/2)) = mat_out;
            kitMedians(ceil(b/2)) = median(median_mat_out, 'all', 'omitnan');
        else
            batMat(:, :, ceil(b/2)) = mat_out;
            batMedians(ceil(b/2)) = median(median_mat_out, 'all', 'omitnan');
        end


    end

    % get average across
    idx = b/2 + 1;
    d.ISC.kitchen_RDM.(frqBand).blockTfa(idx).RDM = mean(kitMat, 3);
    d.ISC.kitchen_RDM.(frqBand).blockTfa(idx).color = [0, 0, 0];
    d.ISC.kitchen_RDM.(frqBand).blockTfa(idx).name = 'average';
    d.ISC.kitchen_RDM.(frqBand).blockTfa(idx).medianISC = mean(kitMedians);

    d.ISC.bathroom_RDM.(frqBand).blockTfa(idx).RDM = mean(batMat, 3);
    d.ISC.bathroom_RDM.(frqBand).blockTfa(idx).color = [0, 0, 0];
    d.ISC.bathroom_RDM.(frqBand).blockTfa(idx).name = 'average';
    d.ISC.bathroom_RDM.(frqBand).blockTfa(idx).medianISC = mean(batMedians);

end
end