function Behaviour_analysis(path_to_data)
% Get the path to the behavioural data, calculates the relevant variables
% and plots figures.
% 
% Kakaei, Schlecht, Hauser 2025
%
% See also Pupil_analysis

if nargin<1
    path_to_data = pwd;
end

for exper = 1:2 % experiment counter
    if exper==1
        arms = 3; % number of bandits
        load(fullfile(path_to_data,'Behav_exp1_dat_anonym'))
    elseif exper==2
        arms = 2; % number of bandits
         load(fullfile(path_to_data,'Behav_exp2_dat_anonym'))
    end

    subs = unique([Behav_dat_anonym.ID]); % subjects' IDs
    %% prepare data
    inotnan = ~isnan(Behav_dat_anonym.Certainty) & ~isnan(Behav_dat_anonym.Value); % has both value and certainty levels
    isfirstchoice = Behav_dat_anonym.choice_index==1; % first response
    Ch1_data_notnan = Behav_dat_anonym(inotnan & isfirstchoice,:); % main data for the analysis

    % mean observed values
    mu1 = Ch1_data_notnan.mu1;
    mu2 = Ch1_data_notnan.mu2;

    chosen_mu = nan(size(mu1)); % mean observed values of the chocen bandit
    notchosen_mu = nan(size(mu1)); % mean observed values of the alternative bandit

    chosen_mu(Ch1_data_notnan.bandit==1) = mu1(Ch1_data_notnan.bandit==1);
    chosen_mu(Ch1_data_notnan.bandit==2) = mu2(Ch1_data_notnan.bandit==2);

    notchosen_mu(Ch1_data_notnan.bandit==1) = mu2(Ch1_data_notnan.bandit==1);
    notchosen_mu(Ch1_data_notnan.bandit==2) = mu1(Ch1_data_notnan.bandit==2);

    dm = chosen_mu-notchosen_mu; % observed value differnece (chosen-not_chosen)
    Ch1_data_notnan.dm = dm;
    %% peaks of P(dm|certainty)
    binsf = linspace(-50,50,200); % Bins of the probability distribution
    P_dm_given_certainty = nan(2,length(binsf),numel(subs)); % P(dm|certainty)
    Certainty = [-1 1]; % uncertain and certain bandits
    peakP = nan(2,numel(subs)); % position of maximum P(dm|certainty)
    for ns = 1:numel(subs)
        Sub_ID = subs(ns); % current subject
        Sub_Dat = Ch1_data_notnan(ismember(Ch1_data_notnan.ID,Sub_ID),:); % subject data

        CC = Sub_Dat.Certainty; % certainty level
        dm = Sub_Dat.dm;
        for ni = 1:2

            [pdensity,xi] = ksdensity(dm(CC==Certainty(ni)),binsf,'Function','pdf'); % prob. density
            [~,idx] = max(pdensity); % position of maximum P(dm|certainty)
            peakP(ni,ns) = xi(idx);
            dp = diff(xi(1:2)); % bin size
            P_dm_given_certainty(ni,:,ns) = pdensity.*dp; % probability density * dp -> P(dm|certainty)
        end
    end
    % consistency of the shift between uncertain and certain bandits over all subjects
    [~,pval,~,stats]  = ttest2(peakP(1,:),peakP(2,:));

    %% figure 1.B (conditional probs)

    f100 = figure(100);
    subplot(1,2,exper)
    MP = mean(P_dm_given_certainty,3,'omitmissing'); % average probability distribution over all subjects
    SP = std(P_dm_given_certainty,0,3,'omitmissing')./sqrt(size(P_dm_given_certainty,3)); % SEM of the distribution over all subjects

    clrs = [0 0.2 0.9;0.7 0 0;0 0.7 0]; % plot colours
    hold on
    plot(binsf,MP(1,:),'Color',clrs(1,:))
    plot(binsf,MP(2,:),'Color',clrs(2,:))
    shade_fig(binsf,MP(1,:)-SP(1,:),MP(1,:)+SP(1,:),clrs(1,:),0.3)
    shade_fig(binsf,MP(2,:)-SP(2,:),MP(2,:)+SP(2,:),clrs(2,:),0.3)

    my_axis_style
    xlabel('\Delta \mu')
    ylabel('P(\Delta \mu | certainty)')
    lg = legend('Uncertain','Certain','Location','northeast');
    lg.Box = 'off';
    lg.AutoUpdate = 'off';
    ylim([0 0.025]);
    yl = ylim;
    xticks(-40:20:40)
    text(0,yl(2),['pval-peaks = ' num2str(pval)])

    %% Figure 1.C (mixed-effect model)
    tmpT = Ch1_data_notnan;
    model =  'isleft~Certainty*dm_left+(Certainty*dm_left|ID)'; %  (logit(P_left) =~β_0  +β_C  C +β_Δμ  Δμ+β_(C⋅Δμ)  (C*Δμ)+ (C*Δμ|subject)+ϵ)

    % left bandit as reference
    isleft = (tmpT.Position==1); % is the selected bandit on the left
    tmpT.isleft = isleft;
    dm = tmpT.dm; % chosen-notchosen
    dm_left = nan(size(dm));
    dm_left(isleft) = dm(isleft); % left - (right or middle) if chosen = left, dm_left = dm
    dm_left(~isleft) = -dm(~isleft); % if notchosen = left, dm_left = -dm

    tmpT.dm_left = dm_left;
  
    lme = fitglme(tmpT,model,'Distribution','binomial') % liner mixed-effect model

    est = lme.Coefficients.Estimate; % betas
    pval = lme.Coefficients.pValue; % pvalues
    CI1 = est-lme.Coefficients.Lower; % lower bound of the confidence interval
    CI2 = lme.Coefficients.Upper-est; % higher bound of the confidence interval

    f101 = figure(101);

    hold on
    bar((exper-1)*2,est(2),'FaceColor',clrs(1,:),'EdgeColor','k','BarWidth',0.3) % beta of Info
    bar((exper-1)*2+0.5,est(3),'FaceColor',clrs(2,:),'EdgeColor','k','BarWidth',0.3) % beta of value
    bar((exper-1)*2+1,est(4),'FaceColor',0.5*ones(1,3),'EdgeColor','k','BarWidth',0.3) % beta of value

    errorbar((exper-1)*2,est(2),CI1(2),CI2(2),'LineWidth',2,'LineStyle','none','Color','k')
    errorbar((exper-1)*2+0.5,est(3),CI1(3),CI2(3),'LineWidth',2,'LineStyle','none','Color','k')
    errorbar((exper-1)*2+1,est(4),CI1(4),CI2(4),'LineWidth',2,'LineStyle','none','Color','k')

    xticks([0.5 2.5])
    xticklabels({'Exp.1','Exp.2'})
    lg = legend('Info','Value','Info*Value','C.I.','Location','northeastoutside');
    lg.Box = 'off';
    lg.AutoUpdate = 'off';
    my_axis_style
    ylabel('\beta')

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