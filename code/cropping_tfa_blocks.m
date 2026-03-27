tf_all2 = tf_all;
durations =[];
Fs = 1

for blk = 1:height(timelock.trial)

    % Get current tf
    blkTf = tf_all2{blk};

    % Get events of current block
    blkTrials = behav(behav.block == blk, :);
    blkEvents = events(blkTrials.trial + (blk - 1)*height(blkTrials));

    % Mark time windows that are good
    taskFreeTPs = false(1, length(blkTf.time));

    for i = 1:height(blkTrials)

        % Define start/end in milliseconds
        t_start = blkTrials.trialOnset(i) - blkTrials.trialOnset(1);  % + 500 ms
        t_end   = t_start + 0.25 + 1;

        durations = [durations, (t_end-t_start)];

        % Convert to samples
        s_start = t_start / Fs;
        s_end   = t_end   / Fs;

        % Mark time points in this window
        windowTPs = blkTf.time < s_end & blkTf.time > s_start;
        taskFreeTPs(windowTPs) = true;

    end

    % Filter for current block
    blkTf.time = blkTf.time(logical(taskFreeTPs));
    blkTf.trial = blk;
    blkTf.powspctrm = blkTf.powspctrm(:, :, :, logical(taskFreeTPs));
    tf_all2{blk} = blkTf;

end

allTPs = nan(20, 4000);
for tfff = 1:length(tf_all2)
    tf_tf = tf_all2{tfff};
    allTPs(tfff, 1:length(tf_tf.time)) = tf_tf.time;
end
figure; imagesc(allTPs);