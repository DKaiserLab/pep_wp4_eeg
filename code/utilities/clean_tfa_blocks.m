function allTfOut = clean_tfa_blocks(cfg)

if ~isfield(cfg, 'threshold'); cfg.threshold = 0.1; end
if ~isfield(cfg, 'nBlks'); cfg.nBlks = 20; end

%% get subjects data
allTf = {[]};
allTfOut = {[]};
for s = 1:length(cfg.subNums)

    % get file name and path
    sub = cfg.subNums(s);
    tfa_path = fullfile(cfg.outputPath, ['sub-', num2str(sub)], 'eeg');
    tfa_filename_all = ['PEP_WP4_EEG', num2str(sub), '_tfa_blocks_all.mat'];

    % load
    allTf{s} = load(fullfile(tfa_path, tfa_filename_all));
end


%% segment and crop chunks (time between task trials) to same length across subjects

for b = 1:cfg.nBlks

    % segment all subjects blocks into chunks
    chunks = cell(cfg.n,1);
    nChunks = zeros(cfg.n,1);

    for s = 1:cfg.n
        blk = allTf{s}.tf_all{b}.time;

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
        allTfOut{s}.tf_all{b}.powspctrm = allTf{s}.tf_all{b}.powspctrm(:, :, :, ...
            ismember(newBlkTimes, allTf{s}.tf_all{b}.time));
        allTfOut{s}.tf_all{b}.time = newBlkTimes;

    end

end


end