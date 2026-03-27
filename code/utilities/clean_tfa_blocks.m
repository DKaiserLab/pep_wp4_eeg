function allTfOut = clean_tfa_blocks(cfg)

if ~isfield(cfg, 'threshold'); cfg.threshold = 0.1; end
if ~isfield(cfg, 'nBlks'); cfg.nBlks = 20; end
if ~isfield(cfg, 'averageChannels'); cfg.averageChannels = true; end
if ~isfield(cfg, 'bands')
    cfg.bands = struct();
    cfg.bands.delta = [1, 3];
    cfg.bands.theta = [4, 7];
    cfg.bands.alpha = [8, 12];
    cfg.bands.beta  = [13, 30];
    cfg.bands.gamma = [31, 70];
end
bandNames = fieldnames(cfg.bands);

%% get subjects data and average frequencies in bands
allTf = {[]};
allTfOut = {[]};
for s = 1:length(cfg.subNums)

    % get file name and path
    sub = cfg.subNums(s);
    tfa_path = fullfile(cfg.outputPath, ['sub-', num2str(sub)], 'eeg');
    tfa_filename_all = ['PEP_WP4_EEG', num2str(sub), '_tfa_blocks_all.mat'];

    % load current tf data
    load(fullfile(tfa_path, tfa_filename_all));

    % loop through blocks
    for b = 1:cfg.nBlks

        pow = tf_all{b}.powspctrm;   % [rpt × chan × freq × time]
        freq = tf_all{b}.freq;

        band_pow = zeros(size(pow,1), size(pow,2), length(bandNames), size(pow,4));

        for bi = 1:length(bandNames)

            band = cfg.bands.(bandNames{bi});

            % find frequency indices
            f_idx = freq >= band(1) & freq <= band(2);

            % average across frequency dimension
            band_pow(:,:,bi,:) = mean(pow(:,:,f_idx,:), 3, 'omitnan');
        end

        % average channels if desired
        if cfg.averageChannels
            band_pow = mean(band_pow, 2, 'omitnan');
        end 

        % store result
        allTf{s}.tf_band{b} = tf_all{b};
        allTf{s}.tf_band{b}.powspctrm = band_pow;       % replace with averaged data
        allTf{s}.tf_band{b}.freq = 1:length(bandNames); % replace freq with band index
        allTf{s}.tf_band{b}.bandlabel = bandNames;      % store names
    end
end


%% segment and crop chunks (time between task trials) to same length across subjects

for b = 1:cfg.nBlks

    % segment all subjects blocks into chunks
    chunks = cell(cfg.n,1);
    nChunks = zeros(cfg.n,1);

    for s = 1:cfg.n
        blk = allTf{s}.tf_band{b}.time;

        d = abs(diff(blk));
        jumpIdx = find(d > cfg.threshold);

        startIdx = [1, jumpIdx + 1];
        endIdx   = [jumpIdx, length(blk)];

        nChunks(s) = length(startIdx);

        % store chunks
        chunks{s} = cell(nChunks(s),1);
        for c = 1:nChunks(s)
            chunks{s}{c} = blk(startIdx(c):endIdx(c));
        end
    end

    % check consistency
    if length(unique(nChunks)) ~= 1
        error('Mismatch in number of chunks across subjects in timeseries %d', b);
    end

    nC = nChunks(1);

    % find minimum chunk length across subjects
    minLengths = zeros(nC,1);

    for c = 1:nC
        lengths = zeros(cfg.n,1);
        for s = 1:cfg.n
            lengths(s) = length(chunks{s}{c});
        end
        minLengths(c) = min(lengths);
    end

    % crop chunks and reconstruct block timeseries
    for s = 1:cfg.n

        newBlkTimes = [];

        for c = 1:nC
            block = chunks{s}{c};

            % crop to minimum length
            croppedBlock = block(1:minLengths(c));

            % concatenate
            newBlkTimes = [newBlkTimes, croppedBlock];
        end


        % store the clean timeseries in the output cell
        %allTfOut{s}.tf_all{b} = ;
        allTfOut{s}.tf_all{b}.powspctrm = allTf{s}.tf_band{b}.powspctrm(:, :, :, ...
            ismember(newBlkTimes, allTf{s}.tf_band{b}.time));
        allTfOut{s}.tf_all{b}.time = newBlkTimes;
        allTfOut{s}.tf_all{b}.bandlabel = allTf{s}.tf_band{b}.bandlabel;
        allTfOut{s}.tf_all{b}.dimord = allTf{s}.tf_band{b}.dimord;

    end

end

% save results
outputFolder = fullfile(cfg.outputPath, 'group_level', 'tfa');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
save(fullfile(outputFolder, 'averageTfAll'), 'allTfOut', '-v7.3')

end