function Pupil_analysis(path_to_data)
% Get the path to the pupil and behaviour data, calculates the relevant
% variables and plots figures.
%
% Kakaei, Schlecht, Hauser 2025
%
% See also Behaviour_analysis

if nargin<1
    path_to_data = pwd;
end

for exper = 1:2 % experiment counter
    if exper==1
        arms = 3; % number of bandits
        load(fullfile(path_to_data,'Pup_exp1_dat_anonym'))
        load(fullfile(path_to_data,'Behav_exp1_dat_anonym'))
    elseif exper==2
        arms = 2; % number of bandits
        load(fullfile(path_to_data,'Pup_exp2_dat_anonym'))
        load(fullfile(path_to_data,'Behav_exp2_dat_anonym'))
    end

    subs = unique([Behav_dat_anonym.ID]); % subjects' IDs
    %% prepare behaviour data

    isfirstchoice = Behav_dat_anonym.choice_index==1; % first response
    Ch1_data = Behav_dat_anonym(isfirstchoice,:); % first choice only

    ihc = Ch1_data.Certainty==1; % high certainty
    ilc = Ch1_data.Certainty==-1; % low certainty

    ihv = Ch1_data.Value==1; % high value
    ilv = Ch1_data.Value==-1; % low value

    %% prepare pupil data
    Alldata = cell2mat(Pup_dat_anonym.Pupil'); % pupil data of all subjects and all trials (timeseries x trials)
    %% Figure 2
    model = 'pup~Certainty+ (1 | ID)'; % linear mixed-effect model
    modeled_indices = 1:100:size(Alldata,1); % model every 100th point (0.1 seconds)

    pval_model = nan(numel(modeled_indices),1); % significance of Beta_Certainty

    dat = Alldata(:,ihc| ilc); % pupil data of only high or low certainty
    tmpT = Ch1_data(ihc| ilc,:); % behaviour data of only high or low certainty

    disp(repmat('|',1,numel(modeled_indices)))
    for nt = 1:numel(modeled_indices) % over time-samples
        fprintf('|')
        iit = modeled_indices(nt); % current time-point
        tmpT.pup = dat(iit,:)'; % pupil data of the current time-point

        lme = fitlme(tmpT,model); % fit the linear model
        pval_model(nt) = lme.Coefficients.pValue(2); % significance of Beta_Certainty
    end
    fprintf('\n')
    % plot
    clrs = [0 0.2 0.9;0.7 0 0]; % figure colours

    f200 = figure(200);
    subplot(1,2,exper)
    ncert = [-1 ,1]; % certainty levels
    hold on
    tab = tmpT;
    t = -3:1e-3:10; t(end) = [];

    for ni = 1:numel(ncert)
        idx = tab.Certainty==ncert(ni); % trials with response curent-bandit
        current_dat = dat(:,idx); % pupil-data of the current bandit
        m = mean(current_dat,2,'omitmissing'); % mean resp
        s = std(current_dat,0,2,'omitmissing')./sqrt(sum(idx)); % SEM
        plot(t,m,'LineWidth',2,'Color',clrs(ni,:)) % plot mean pupil size over time
        shade_fig(t,m-s,m+s,clrs(ni,:),0.3); % plot SEM
    end
    lg = legend('Uncertain','SEM','Certain','Location','northeast');
    lg.Box = 'off';
    lg.AutoUpdate = 'off';

    ylim([-1 1])
    tmpy = ylim;

    % show the significant p-values
    ip = pval_model<0.05; % significant
    idx = modeled_indices(ip); % indices of the significant time points
    tmpy = ones(size(t))*tmpy(2); % where to put the asterisks

    plot(t(idx),tmpy(idx)-0.1,'*','Color','k')
    xlabel('trial onset (s)')
    ylabel('mean pupil size (Z-scored)')
    title(['Experiment ' num2str(exper)])
    my_axis_style
    %% Figure2 - subplot Time window [-1 0]
    TOI = [-1 0]; % time window of interest
    tmpT = Ch1_data(ihc | ilc,:); % only with certainty disparity

    tt = -3:1e-3:10; tt(end) = []; % time vector

    iit1 = find(tt<TOI(1),1,'last'); % first time point of interest
    iit2 = find(tt<TOI(2),1,'last'); % last time point of interest
    MPup = mean(Alldata(iit1:iit2,ihc | ilc),1)'; % mean_pupil size within TOI

    Cert = (tmpT.Certainty); % certainty levels

    MP_TOI = nan(2,numel(subs));
    for ns = 1:numel(subs)
        tmp1 = MPup(Cert==-1 & ismember(tmpT.ID,subs(ns)),:);
        tmp2 = MPup(Cert==1 & ismember(tmpT.ID,subs(ns)),:);

        MP_TOI(1,ns) = mean(tmp1,"omitmissing"); % average over trials
        MP_TOI(2,ns) = mean(tmp2,"omitmissing"); % average over trials
    end

    f201 = figure(201);
    subplot(1,2,exper)
    hold on
    b = bar(-0.22,mean(MP_TOI(1,:),2),0.4); % average over subjects
    b.FaceColor = clrs(1,:);
    b.FaceAlpha = 1;
    b.EdgeColor = 'k';
    errorbar(-0.22,mean(MP_TOI(1,:),2),std(MP_TOI(1,:),0,2)/sqrt(size(MP_TOI,2)),'LineWidth',2,'Color','k') % SEM over subjects

    b = bar(0.22,mean(MP_TOI(2,:),2),0.4);
    b.FaceColor = clrs(2,:);
    b.FaceAlpha = 1;
    b.EdgeColor = 'k';
    errorbar(0.22,mean(MP_TOI(2,:),2),std(MP_TOI(2,:),0,2)/sqrt(size(MP_TOI,2)),'LineWidth',2,'Color','k')
    lg = legend('Uncertain','SEM','Certain','Location','northeast');
    lg.Box = 'off';
    lg.AutoUpdate = 'off';
    ylabel('mean pupil size (Z-scored)')
    xticks([]);
    my_axis_style
    title(['Experiment ' num2str(exper)])
    ylim([-0.7 -0.2])
    yticks(-0.6:0.2:-0.2)

    % linear mixed-model
    modeled_indicestmp = iit1:100:iit2; % model every 100th point
    w = length(modeled_indicestmp); % window size

    % certainty based, controlled for autocorr
    model = 'pup~Certainty+(1|ID)+(1|time)';
    tmpT = Ch1_data(ihc | ilc,:); % behaviour data only with certainty disparity
    dat = Alldata   (:,ihc | ilc); % pupil data only with certainty disparity
    tmpT.pup = nan(size(tmpT,1),1);
    tmpT.time = nan(size(tmpT,1),1);
    longTcert = repmat(tmpT,w,1); % initialise (will be overwritten)

    ii = 1;
    for ntr = 1:size(tmpT,1) % over trials
        iit = modeled_indicestmp;
        longTcert(ii:(ii+w-1),:) = repmat(tmpT(ntr,:),w,1); % every parameter is repeated except the pupil data
        longTcert.pup(ii:(ii+w-1)) = dat(iit,ntr); % aggrigate data of the window
        longTcert.time(ii:(ii+w-1)) = iit; % time points
        ii = ii+w;
    end
    disp(['Figure 2 - Experiment ' num2str(exper)])
    lme = fitlme(longTcert,model) % fit the linear model

    pval_model_windowed_info = lme.Coefficients.pValue(2); % significance of Beta_Certainty
    %% Figure 3
    model = 'pup~Value+ (1 | ID)'; % linear mixed-effect model
    modeled_indices = 1:100:size(Alldata,1); % model every 100th point (0.1 seconds)

    pval_model = nan(numel(modeled_indices),1); % significance of Beta_Certainty

    dat = Alldata(:,ihv| ilv); % pupil data of only high or low Value
    tmpT = Ch1_data(ihv| ilv,:); % behaviour data of only high or low Value

    disp(repmat('|',1,numel(modeled_indices)))
    for nt = 1:numel(modeled_indices) % over time-samples
        fprintf('|')
        iit = modeled_indices(nt); % current time-point
        tmpT.pup = dat(iit,:)'; % pupil data of the current time-point

        lme = fitlme(tmpT,model); % fit the linear model
        pval_model(nt) = lme.Coefficients.pValue(2); % significance of Beta_Certainty
    end
    fprintf('\n')
    % plot
    clrs = [0 0.2 0.9;0.7 0 0]; % figure colours

    f300 = figure(300);
    subplot(1,2,exper)
    nval = [-1 ,1]; % Value levels
    hold on
    tab = tmpT;
    t = -3:1e-3:10; t(end) = [];

    for ni = 1:numel(nval)
        idx = tab.Value==nval(ni); % trials with response curent-bandit
        current_dat = dat(:,idx); % pupil-data of the current bandit
        m = mean(current_dat,2,'omitmissing'); % mean resp
        s = std(current_dat,0,2,'omitmissing')./sqrt(sum(idx)); % SEM
        plot(t,m,'LineWidth',2,'Color',clrs(ni,:)) % plot mean pupil size over time
        shade_fig(t,m-s,m+s,clrs(ni,:),0.3); % plot SEM
    end
    lg = legend('High-Value','SEM','Low-Value','Location','northeast');
    lg.Box = 'off';
    lg.AutoUpdate = 'off';

    ylim([-1 1])
    tmpy = ylim;

    % show the significant p-values
    ip = pval_model<0.05; % significant
    idx = modeled_indices(ip); % indices of the significant time points
    tmpy = ones(size(t))*tmpy(2); % where to put the asterisks

    plot(t(idx),tmpy(idx)-0.1,'*','Color','k')
    xlabel('trial onset (s)')
    ylabel('mean pupil size (Z-scored)')
    title(['Experiment ' num2str(exper)])
    my_axis_style
     %% Figure3 - subplot Time window [-1 0]
    TOI = [-1 0]; % time window of interest
    tmpT = Ch1_data(ihv| ilv,:); % only with certainty disparity

    tt = -3:1e-3:10; tt(end) = []; % time vector

    iit1 = find(tt<TOI(1),1,'last'); % first time point of interest
    iit2 = find(tt<TOI(2),1,'last'); % last time point of interest
    MPup = mean(Alldata(iit1:iit2,ihv | ilv),1)'; % mean_pupil size within TOI

    VAL = (tmpT.Value); % certainty levels

    MP_TOI = nan(2,numel(subs));
    for ns = 1:numel(subs)
        tmp1 = MPup(VAL==-1 & ismember(tmpT.ID,subs(ns)),:);
        tmp2 = MPup(VAL==1 & ismember(tmpT.ID,subs(ns)),:);

        MP_TOI(1,ns) = mean(tmp1,"omitmissing"); % average over trials
        MP_TOI(2,ns) = mean(tmp2,"omitmissing"); % average over trials
    end

    f301 = figure(301);
    subplot(1,2,exper)
    hold on
    b = bar(-0.22,mean(MP_TOI(1,:),2),0.4); % average over subjects
    b.FaceColor = clrs(1,:);
    b.FaceAlpha = 1;
    b.EdgeColor = 'k';
    errorbar(-0.22,mean(MP_TOI(1,:),2),std(MP_TOI(1,:),0,2)/sqrt(size(MP_TOI,2)),'LineWidth',2,'Color','k') % SEM over subjects

    b = bar(0.22,mean(MP_TOI(2,:),2),0.4);
    b.FaceColor = clrs(2,:);
    b.FaceAlpha = 1;
    b.EdgeColor = 'k';
    errorbar(0.22,mean(MP_TOI(2,:),2),std(MP_TOI(2,:),0,2)/sqrt(size(MP_TOI,2)),'LineWidth',2,'Color','k')
    lg = legend('Low value','SEM','High value','Location','northeast');
    lg.Box = 'off';
    lg.AutoUpdate = 'off';
    ylabel('mean pupil size (Z-scored)')
    xticks([]);
    my_axis_style
    title(['Experiment ' num2str(exper)])
    ylim([-0.7 -0.2])
    yticks(-0.6:0.2:-0.2)

    % linear mixed-model
    modeled_indicestmp = iit1:100:iit2; % model every 100th point
    w = length(modeled_indicestmp); % window size

    % certainty based, controlled for autocorr
    model = 'pup~Value+(1|ID)+(1|time)';
    tmpT = Ch1_data(ihv | ilv,:); % behaviour data only with certainty disparity
    dat = Alldata   (:,ihv | ilv); % pupil data only with certainty disparity
    tmpT.pup = nan(size(tmpT,1),1);
    tmpT.time = nan(size(tmpT,1),1);
    longTcert = repmat(tmpT,w,1); % initialise (will be overwritten)

    ii = 1;
    for ntr = 1:size(tmpT,1) % over trials
        iit = modeled_indicestmp;
        longTcert(ii:(ii+w-1),:) = repmat(tmpT(ntr,:),w,1); % every parameter is repeated except the pupil data
        longTcert.pup(ii:(ii+w-1)) = dat(iit,ntr); % aggrigate data of the window
        longTcert.time(ii:(ii+w-1)) = iit; % time points
        ii = ii+w;
    end
    disp(['Figure 3 - Experiment ' num2str(exper)])
    lme = fitlme(longTcert,model) % fit the linear model

    pval_model_windowed_info = lme.Coefficients.pValue(2); % significance of Beta_Certainty
    %%
end
end
%% functions
function my_axis_style
% adjust fonts and line width
%
% Ehsan Kakaei 2021

fs = 15;
axis square
box off
ax = gca;
ax.LineWidth = 2;
ax.FontSize = fs;
ax.FontName = 'Arial';

end
function shade_fig(x,ylow,yhigh,c,a)
% creates a shaded area on current figure
%
% SHADE_FIG(X,YLOW,YHIGH,C,A) creates a shaded area with higher limit YHIGH
% and lower limit YLOW on x-coordinates X, with the color C and
% transparency A.
%
% Ehsan Kakaei 2021

x = x(:);
ylow = ylow(:);
yhigh = yhigh(:);
f = fill([x;flipud(x)],[ylow;flipud(yhigh)],c);
f.FaceAlpha = a;
f.EdgeColor = 'none';
end