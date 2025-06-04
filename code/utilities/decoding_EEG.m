%% Housekeeping
clear 
clc
close all

%% Define Parameters
ss=[102 103 104 106 109 110 111 112 113 115 116 117 120 122];
sn=length(ss);
s_idx=1;

output_dir=fullfile(pwd, '..', '..', 'derivatives');
smoothing_window=6;

%cd('C:\MATLAB\Toolboxes\CoSMoMVPA-master\mvpa')%'C:\temp\Scene
%Clips 2\Toolboxes\CoSMoMVPA-master\mvpa\'
%addoath 'C:\MATLAB\Toolboxes\fieldtrip-20220104\
ft_defaults;

%% Decoding
for s=ss% for all subjects
    
    filepath = fullfile(pwd, '..', '..', 'derivatives', ['sub-', num2str(s)], 'eeg', ['PEP_WP4_EEG', num2str(s), '_timelock', '.mat']);
    load(filepath);
    
    % get time info
    res.time=timelock.time;

    % convert to cosmo
    ds=cosmo_meeg_dataset(timelock);
    clear timelock

    % add targets
    ds.sa.targets=ds.sa.trialinfo(:,1);

    % add chunks
    nch=100; % für pairwise:20

    ds.sa.chunks=[1:length(ds.sa.targets)]';
    ds.sa.chunks=cosmo_chunkize(ds,nch);
    
    
    % conduct decoding analysis at every time point
    
    for tp=1:max(ds.fa.time)% for all time points

        display(['Subject ' num2str(s_idx) ' of ' num2str(sn) '. Time point ' num2str(tp) ' of ' num2str(max(ds.fa.time)) '.']);

        ds_tp=cosmo_slice(ds,ismember(ds.fa.time,tp),2);

        % do pca

        n_feat=length(unique(ds_tp.fa.chan));
        [coeff,x,LATENT,~,x_exp,mu]=pca(ds_tp.samples);
        for ccx=1:length(x_exp)
           if sum(x_exp(1:ccx))>=99
              n_feat=ccx;
             break
          end
        end
        ds_tp.samples=x(:,1:n_feat);
        ds_tp.fa.chan=1:n_feat;
        ds_tp.fa.time=repmat(tp,1,n_feat);
        ds_tp.a.fdim.values{1}=1:n_feat;

        ds_class=ds_tp;

        % assign targets (i.e. scene trigger codes) to categories
        % 1-50 bathrooms
        % 51-100 kitchen

        bathroom_logical=ismember(ds.sa.targets,[1:50]);
        kitchen_logical=ismember(ds.sa.targets,[51:100]);

        ds_class.sa.targets(bathroom_logical,1)=1;
        ds_class.sa.targets(kitchen_logical,1)=2;

        % decoding settings
        args.classifier=@cosmo_classify_lda;
        args.partitions=cosmo_nchoosek_partitioner(ds_class,1);

        % run decoding
        acc=cosmo_crossvalidation_measure(ds_class,args);
        res.dec_acc(s_idx,tp)=acc.samples;
        res.order(s_idx) = s;

    end % time points

s_idx=s_idx+1;
end % subjects

% smooth the decoding curve with a rolling average
res.dec_acc=smoothdata(res.dec_acc,2,'movmean',smoothing_window);

file_name=['PEP_WP4_EEG_decoding_accuracy',num2str(smoothing_window),'_TP_RA.mat'];

% save decoding data to file
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

cd(output_dir);
save(file_name,'res');

%% plot decoding results

dec_acc_mean=squeeze(mean(res.dec_acc,1));

figure();
decoding_plot=plot(res.time, dec_acc_mean(:));
chance=1/2;
     
yline(chance,'k--');
xline(0,'k');
title(['pep_wp4_eeg decoding']);

