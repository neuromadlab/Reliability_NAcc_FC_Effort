% ======================================================================= %
% How reliable is the link between effort and nucleus accumbens
% connectivity at rest? A psychometric evaluation of repeated assessments
%
% Functional connectome fingerprinting (Results; Figure 3B, 3C)
%
% Input:  data/fingerprinting_similarity_matrix.csv
%         Pearson r between the gray-matter NAcc FC maps of every participant
%         in S6 (rows) and S7 (columns); participants in the same order.
% Output: output/fingerprinting_results.mat, output/Figure3B/3C_*.png
%
% Run the whole file (e.g. run('code/fingerprinting.m')). Requires the Statistics and Machine Learning
% Toolbox. Figure 3B uses the 'plasma' colormap if it is on the path
% (otherwise parula).
% ======================================================================= %

clearvars; clc; close all;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    repo_dir = pwd;
else
    repo_dir = fileparts(script_dir);
end
data_dir   = fullfile(repo_dir,'data');
output_dir = fullfile(repo_dir,'output');

sim_file = fullfile(data_dir,'fingerprinting_similarity_matrix.csv');
if ~isfile(sim_file)
    error(['Cannot find %s\nRun the script from the extracted repository ' ...
           '(not from inside the zip), or cd to the repository root.'],sim_file);
end
if ~exist(output_dir,'dir'), mkdir(output_dir); end


%% Similarity matrix ------------------------------------------------------

opts = detectImportOptions(sim_file,'VariableNamingRule','preserve');
opts.VariableNamesLine = 1;
opts.DataLines = [2 Inf];
opts = setvartype(opts,1,'string');
opts = setvartype(opts,2:numel(opts.VariableNames),'double');
T = readtable(sim_file,opts);

if ~strcmp(T.Properties.VariableNames{1},'ID')
    error('First column must be ID (anonymised data/ file, not the internal one).');
end
IDs = T{:,1};
corr_matrix = T{:,2:end};
N = numel(IDs);

if ~isequal(size(corr_matrix),[N N]) || ...
        ~isequal(string(T.Properties.VariableNames(2:end))',IDs)
    error('Rows (S6) and columns (S7) must list the same participants in the same order.');
end


%% Identification accuracy and similarity ---------------------------------

[~,max_idx] = max(corr_matrix,[],2);
accuracy = mean((1:N)' == max_idx);

within_subj_r  = diag(corr_matrix);
between_subj_r = corr_matrix(~eye(N));

corr_clipped = corr_matrix;
corr_clipped(corr_clipped >= 1)  = 0.999999;
corr_clipped(corr_clipped <= -1) = -0.999999;
z_corr_matrix  = atanh(corr_clipped);
within_subj_z  = diag(z_corr_matrix);
between_subj_z = z_corr_matrix(~eye(N));

mean_r_within  = mean(within_subj_r,'omitnan');
mean_r_between = mean(between_subj_r,'omitnan');
mean_z_within  = mean(within_subj_z,'omitnan');
mean_z_between = mean(between_subj_z,'omitnan');

% Permutation test: shuffle the S7 participants
rng(2026);
nPerm = 1000;
perm_acc = zeros(nPerm,1);
for p = 1:nPerm
    [~,max_idx_perm] = max(corr_matrix(:,randperm(N)),[],2);
    perm_acc(p) = mean((1:N)' == max_idx_perm);
end
p_perm = (sum(perm_acc >= accuracy) + 1) / (nPerm + 1);

chance_level = 1/N;
n_successes  = round(accuracy*N);
p_binomial   = 1 - binocdf(n_successes-1,N,chance_level);

fprintf('\nParticipants: %d\n',N);
fprintf('Fingerprinting accuracy = %.2f%% (chance = %.2f%%)\n',accuracy*100,chance_level*100);
fprintf('Mean within / between Pearson r = %.3f / %.3f\n',mean_r_within,mean_r_between);
fprintf('Mean within / between Fisher z  = %.3f / %.3f\n',mean_z_within,mean_z_between);
fprintf('Permutation p = %.4f (%d permutations); binomial p = %.4g\n',p_perm,nPerm,p_binomial);

save(fullfile(output_dir,'fingerprinting_results.mat'), ...
    'corr_matrix','z_corr_matrix','accuracy','within_subj_r','between_subj_r', ...
    'within_subj_z','between_subj_z','mean_r_within','mean_r_between', ...
    'mean_z_within','mean_z_between','perm_acc','p_perm','p_binomial','IDs');


%% Figure 3B: participant-by-participant similarity matrix ---------------

fig = figure('Color','w','Position',[100 100 800 800]);
imagesc(corr_matrix,[0 0.5]);
axis square;
if exist('plasma','file') == 2
    colormap(plasma);
else
    colormap(parula);
end
colorbar;
xlabel('Session 7 participants');
ylabel('Session 6 participants');
title(sprintf('Identification accuracy = %.2f%%',accuracy*100));
ticks = unique([1 5:5:N]);
set(gca,'FontSize',16,'XTick',ticks,'YTick',ticks);
saveas(fig,fullfile(output_dir,'Figure3B_fingerprinting_matrix.png'));


%% Figure 3C left: empirical cumulative distributions --------------------

[f_w,x_w,flo_w,fup_w] = ecdf(within_subj_z,'alpha',0.05);
[f_b,x_b,flo_b,fup_b] = ecdf(between_subj_z,'alpha',0.05);

fig = figure('Color','w','Position',[100 100 650 540]);
hold on;
h_b  = plot(x_b,f_b,'r-','LineWidth',2);
h_w  = plot(x_w,f_w,'b-','LineWidth',2);
h_wc = plot(x_w,flo_w,'b--','LineWidth',0.8);
plot(x_w,fup_w,'b--','LineWidth',0.8,'HandleVisibility','off');
h_bc = plot(x_b,flo_b,'r--','LineWidth',0.8);
plot(x_b,fup_b,'r--','LineWidth',0.8,'HandleVisibility','off');
xlabel('Fisher z-transformed correlation');
ylabel('Empirical cumulative distribution');
legend([h_b h_w h_wc h_bc],{'Off-diagonal','Diagonal','95% CI (within)','95% CI (between)'}, ...
    'Location','southeast','Box','off');
grid on; box off;
set(gca,'FontSize',18);
xl = xlim;
text(xl(1)+0.05*diff(xl),0.15,sprintf('%.2f',mean_z_between),'Color','r','FontSize',18,'FontWeight','bold');
text(xl(2)-0.25*diff(xl),0.85,sprintf('%.2f',mean_z_within),'Color','b','FontSize',18,'FontWeight','bold');
saveas(fig,fullfile(output_dir,'Figure3C_ECDF.png'));


%% Figure 3C right: mean within vs between Fisher z (+/- SEM) ------------

sem = @(v) std(v,'omitnan') / sqrt(sum(isfinite(v)));

fig = figure('Color','w','Position',[100 100 520 520]);
b = bar([mean_z_within mean_z_between],'FaceColor','flat');
b.CData = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
hold on;
errorbar(1:2,[mean_z_within mean_z_between],[sem(within_subj_z) sem(between_subj_z)], ...
    'k','LineStyle','none','LineWidth',1.5);
set(gca,'XTick',1:2,'XTickLabel',{'Diagonal','Off-diagonal'},'FontSize',14, ...
    'YLim',[0 0.5],'YTick',0:0.1:0.5);
ylabel('Mean Fisher z-transformed correlation');
box off; grid on;
saveas(fig,fullfile(output_dir,'Figure3C_mean_similarity.png'));