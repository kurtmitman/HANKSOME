clear
clc
close all

% Load solution
if(exist('exp_both_41_lev16.mat') == 0)
load solution_figures  
else
load exp_both_41_lev16
end    
    
mod_slide=1;
col_poss_20=lines(20);
col_poss=[
   col_poss_20(6,:); col_poss_20(1,:);col_poss_20(2,:)];

Lins={'-','-','-'};
picture_scale=0.5;
Tirf=100;

%% Figure 1 

font_s=16;
font_l=24;
figure(1)

axpos = get(gca,'pos');
orient landscape
hold on
plot(0:TT-1,-0.25*shock_b_agg/aggregates_initial.y,'LineWidth',3,'color',col_poss(1,:),'LineStyle','--')
plot(shock_start_t-1:TT,-0.25*[shock_b_agg(shock_start_t-1)/aggregates_initial.y;new_shock_b_agg(1:TT-shock_start_t+1)/aggregates_initial.y],'LineWidth',3,'color',col_poss(2,:))
line([shock_start_t-1 shock_start_t-1],[0 1],'LineWidth',3,'Color',[1 .2 .4],'LineStyle',':')

ylabel('% of annual GDP','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)

set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Foreign Credit Supply}','fontname','times','fontsize',font_l,'Interpreter','Latex')

grid on
l=legend('Credit Expansion','Unexpected Contraction','location','southeast');
set(l,'Interpreter','Latex')
legend boxoff

%% Figure 2
figure(2)
axpos = get(gca,'pos');
 set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-0.02])
 set(gcf,'Position',[50 50 picture_scale*1800 picture_scale*700*.5./.3]);
 set(gcf, 'PaperPositionMode', 'auto');
 orient landscape
 
 
subplot(2,2,1)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([aggregates_initial.c_h; c_t(1:Tirf-1)])-log(aggregates_initial.c_h)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')


axis([0 Tirf-1 0 10]);
ylabel('% dev. from initial ss','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Consumption of Home Goods}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,2)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

cf=params.omega./(1-params.omega)*[tot'].*c_t;
cf_agg_ss2=aggregates_initial.c_h*T_initial*params.omega./(1-params.omega);
plot(0:Tirf-1,100*(log([cf_agg_ss2; cf(1:Tirf-1)])-log(cf_agg_ss2)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})

grid on

line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

axis([0 Tirf-1 -5 15]);
ylabel('% dev. from initial ss','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Consumption of Foreign Goods}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,3)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([aggregates_initial.y; y_agg(1:Tirf-1)])-log(aggregates_initial.y)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{1})

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 -5 10]);
ylabel('%  dev. from initial ss','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Output}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,4)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

inv=k_agg(2:TT)-(1-params.delta)*k_agg(1:TT-1);
plot(0:Tirf-1,100*(([params.delta*aggregates_initial.k; inv(1:Tirf-1)])-params.delta*aggregates_initial.k)./(params.delta*aggregates_initial.k),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 0 20]);
ylabel('% dev. from initial ss','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Investment}','fontname','times','fontsize',font_l,'Interpreter','Latex')

%% Figure 3

% Terms trade real wage Nominal exc rate real interest rate
figure(3)
axpos = get(gca,'pos');
 set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-0.02])
 set(gcf,'Position',[50 50 picture_scale*1800*1. picture_scale*700*1.*.5./.3]);
 set(gcf, 'PaperPositionMode', 'auto');
 orient landscape
 
subplot(2,2,1)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([T_initial; tot(1:Tirf-1)'])-log(T_initial)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{2})

line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

hold on
axis([0 Tirf-1 -5 10]);
grid on
ylabel('% dev. from initial ss','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Terms of Trade}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,2)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([aggregates_initial.w; w_nr_guess(1:Tirf-1)])-log(aggregates_initial.w)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})

grid on

line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

axis([0 Tirf-1 0 3]);

ylabel('% dev. from initial ss','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Real Wage}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,3)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')


plot(0:Tirf-1,100*(log([1/T_initial; tot(1:Tirf-1)'.^(-1)])-log(1/T_initial)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 -10 10]);
ylabel('%  dev. from initial ss','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Nominal Exchange Rate}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,4)
hold on
line([-100 100],100*[params.r_ss params.r_ss],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(([aggregates_initial.r_k-params.delta; r_t(1:Tirf-1)])),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 0 4]);
ylabel('% Level','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Real Interest Rate }','fontname','times','fontsize',font_l,'Interpreter','Latex')


%% Define Foreign Consumption, Investment and Real wage

Tirf=30;

Consumption_foreign_fixed= params.omega./(1-params.omega)*tot_fixed.*c_t_fixed;
Consumption_foreign_exp= params.omega./(1-params.omega)*tot_exp'.*c_t_exp;

k_agg0 = k_agg(shock_start_t-1);
ch_agg0 = c_t_agg(shock_start_t);
cf_agg0 = params.omega./(1-params.omega)*tot(shock_start_t).*c_t(shock_start_t);
y_agg0 = y_agg(shock_start_t);
inv0 = k_agg(shock_start_t)-(1-params.delta)*k_agg(shock_start_t-1);
k_agg_past_exp=[k_agg0;k_agg_exp(1:end-1)];
k_agg_past_fixed=[k_agg0;k_agg_fixed(1:end-1)];
Investment_exp =k_agg_exp-(1-params.delta)*k_agg_past_exp;
Investment_fixed =k_agg_fixed-(1-params.delta)*k_agg_past_fixed;
tot0=tot(shock_start_t);
w0 = w_nr_guess(shock_start_t);
r0 = r_t(shock_start_t+1);

%% Figure 4

figure(4)
axpos = get(gca,'pos');
 set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-0.02])
 set(gcf,'Position',[50 50 picture_scale*1800 picture_scale*700*.5./.3]);
 set(gcf, 'PaperPositionMode', 'auto');
 orient landscape
 
 
subplot(2,2,1)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([ch_agg0; c_t_exp(1:Tirf-1)])-log(ch_agg0)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})
plot(0:Tirf-1,100*(log([ch_agg0; c_t_fixed(1:Tirf-1)])-log(ch_agg0)),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--')

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')


axis([0 Tirf-1 -10 5]);
ylabel('% dev. from transition path','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Consumption of Home Goods}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,2)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')


plot(0:Tirf-1,100*(log([cf_agg0; Consumption_foreign_exp(1:Tirf-1)])-log(cf_agg0)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})
plot(0:Tirf-1,100*(log([cf_agg0; Consumption_foreign_fixed(1:Tirf-1)])-log(cf_agg0)),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--')

grid on

line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')


axis([0 Tirf-1 -25 5]);
ylabel('% dev. from transition path','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Consumption of Foreign Goods}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,3)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([y_agg0; y_agg_exp(1:Tirf-1)])-log(y_agg0)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{1})
plot(0:Tirf-1,100*(log([y_agg0; y_agg_fixed(1:Tirf-1)])-log(y_agg0)),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--')

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 -5 5]);
ylabel('%  dev. from transition path','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Output}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,4)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot1=plot(0:Tirf-1,100*(([inv0; Investment_exp(1:Tirf-1)])-inv0)./(inv0),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk});
plot2=plot(0:Tirf-1,100*(([inv0; Investment_fixed(1:Tirf-1)])-inv0)./(inv0),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--');

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 -20 10]);
ylabel('% dev. from transition path','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Investment}','fontname','times','fontsize',font_l,'Interpreter','Latex')
l=legend([plot1 plot2],'Flexible','Fixed','location','southeast');
set(l,'Interpreter','Latex')
legend boxoff

%% Figure 5
% Terms trade real wage Nominal exc rate real interest rate

figure(5)
axpos = get(gca,'pos');
 set(gca,'pos',[axpos(1) axpos(2) axpos(3) axpos(4)-0.02])
 set(gcf,'Position',[50 50 picture_scale*1800*1. picture_scale*700*1.*.5./.3]);
 set(gcf, 'PaperPositionMode', 'auto');
 orient landscape
 
 
subplot(2,2,1)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([tot0; tot_exp(1:Tirf-1)'])-log(tot0)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{2})
plot(0:Tirf-1,100*(log([tot0; tot_fixed(1:Tirf-1)])-log(tot0)),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--')

line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

hold on
axis([0 Tirf-1 -15 2]);
grid on
ylabel('% dev. from transition path','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Terms of Trade}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,2)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')


plot(0:Tirf-1,100*(log([w0; w_nr_guess_exp(1:Tirf-1)])-log(w0)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})
plot(0:Tirf-1,100*(log([w0; w_nr_guess_fixed(1:Tirf-1)])-log(w0)),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--')

grid on

line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

axis([0 Tirf-1 -4 2]);
ylabel('% dev. from transition path','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Real Wage}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,3)
hold on
line([-100 100],[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot(0:Tirf-1,100*(log([1./tot0; tot_exp(1:Tirf-1)'.^(-1)])-log(tot0)),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk})
plot(0:Tirf-1,100*(log([1; ones(Tirf-1,1)])-log(1)),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--')

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 -2 10]);
ylabel('%  dev. from transition path','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Nominal Exchange Rate}','fontname','times','fontsize',font_l,'Interpreter','Latex')

subplot(2,2,4)
hold on
line([-100 100],100*[0 0],'LineWidth',1,'Color',[1 .2 .4],'LineStyle',':')

plot1=plot(0:Tirf-1,100*(([0; r_t_exp(2:Tirf)-r0])),'LineWidth',3,'color',col_poss(2,:),'LineStyle',Lins{kk});
plot2=plot(0:Tirf-1,100*(([0; r_t_fixed(2:Tirf)-r0])),'LineWidth',3,'color',col_poss(3,:),'LineStyle','--');

grid on
line([1 1],[-100 100],'LineWidth',2,'Color',[.5 .5 .5],'LineStyle',':')

grid on
axis([0 Tirf-1 -1 3]);
ylabel('% Level','fontname','times','fontsize',font_s)
xlabel('Quarters','fontname','times','fontsize',font_s)
set(gca,'box','on','fontname','times','fontsize',font_s)
title('\textbf{Real Interest Rate }','fontname','times','fontsize',font_l,'Interpreter','Latex')
l=legend([plot1 plot2],'Flexible','Fixed','location','southeast');
set(l,'Interpreter','Latex')
legend boxoff