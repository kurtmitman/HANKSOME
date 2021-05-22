%% Two-Good Open Economy, Steady state transition - 31 Jan 19

clear
clc
close all

%% Add to path folders with additional files

addpath([pwd,'/Steady_State_Ayagari'])
addpath([pwd,'/Transition'])
addpath([pwd,'/Exogenous_files'])

%% Target moments

target.gr_inc = 7.9; % This is p-c gr income in thousands of euro, hungary (I1)
target.deb_inc = .603; % Median d-i ratio, F4
target.gdp = 12.65/4; 
target.nw = 26.2; % Net wealth, per capita A3, th. of eur.
target.nfa = -6.5*.69; %NFA per cap., th ou eur, from eurostat (TOT = 65 bn eur)
target.cap = -target.nfa+target.nw;
target.cap_gdp = target.cap./target.gdp;
target.nfa_gdp = target.nfa./target.gdp;
params.cur_acc = 0.074; % Expansion of current account/GDP in the initial period

target.nw_quintiles = [.9 13 26.2 48.6 84.1 160.6];

%% I. Parameter Input

Ind_sol = 1; % Set to unity if you don't want to calibrate beta

params.alpha        = 0.33; % capital-income share
params.gamma        = 1;    % CRRA (must be > 0)
params.theta        = 1/params.gamma; % Restrict intra-EOS to equal 1/ies
params.delta        = 0.08/4;% depreciation rate
params.omega        = .4;	% Share of foreign goods in consumption
params.z_ss         = 1;	% SS-Productivity in T-sector
params.T_ss         = 1;
params.nr           =.96; % Schmidt-Grohe and Uribe Evidence for Argentina
params.varphi       =.5; % Frish Elasticity - Kurt Paper with Marcus and Iouri
params.epsilon_f       =10; % Substitution across varieties
params.epsilon_w       =10; % Substitution across varieties
params.theta_nr       =100; % epsilon/theta=.1 as HANK


% Process for Labor Endowment
params.nx   = 7;        % Number of points for process
nx          = params.nx;
params.rho  = .97; %Persistence of AR(1)
rho         = params.rho;
sigy        = .84*sqrt(1-params.rho^2); %Variance of innovations of AR(1)
ns          = nx;
params.ns   = params.nx;
[params.ex,params.piex] = rouwenhorst(params.nx,0,rho,sigy);
params.ex               = exp(params.ex); % Labor endowment process
temp = params.piex^1000000;
params.exinv =temp(1,:);    % Ergodic distribution of exogenous process
ydist = params.exinv;
params.ex = params.ex/(params.exinv*params.ex); % Normalize mean to 1
params.ex_mean=params.exinv*params.ex; % Mean Labor Endowment
params.l_bar  =1;


%Steady state quarterly K/Y ratio = 10.26, normalize H=1
params.k_y_ss   = target.cap_gdp;  % Capital-Output ratio. Target value
params.l_ss     = params.l_bar*params.ex_mean;   % Normalize tot. hours to 1. 
params.k_ss     = params.l_ss*(params.z_ss*params.k_y_ss)...
    ^(1/(1-params.alpha)); % SS-capital stock
params.y_ss  = params.z_ss*params.k_ss^params.alpha*...
    params.l_ss^(1-params.alpha); % SS T-Output
params.r_ss  =(params.epsilon_f -1)./params.epsilon_f .* params.alpha*params.z_ss*params.k_ss^(params.alpha-1)*...
    params.l_ss^(1-params.alpha) - params.delta; % SS return on capital, net of depreciation
params.pi_ss  = 1/params.epsilon_f*params.y_ss; %% SS profits
params.w_ss =(params.epsilon_f -1)./params.epsilon_f.*(1-params.alpha)*params.y_ss/params.l_ss; % SS wage

%Borrowing constraint [Units of mean quarterly income]
params.bmin     = -10^(-7)*params.k_y_ss*params.y_ss;

% External Asset Supply (Negative if ROW is lending to our country)
params.b_agg_ss = target.nfa_gdp*params.y_ss; 
 
% Consumption of home good
params.c_h =  (1+(params.omega/(1-params.omega))).^(-1)...
            .*(params.r_ss.*(params.b_agg_ss+params.k_ss) +...
            params.l_ss.*params.w_ss + params.pi_ss);
 
% Labor Disutility        
params.chi_dis=(1-params.omega)* params.w_ss.*(params.epsilon_w-1)./...
    (params.epsilon_w.*...
    params.c_h.*params.l_ss.^(params.varphi));


params.Amax = -params.bmin*999/10^(-7); % maximum asset holdings


% 2. technical parameters to solve model
curv        = 3;
if mod(curv,2)==0
    error('curv must be an odd number')
end
linsp_min=-.2;
Initial_value=10+abs(linsp_min)*10+1;

for s=1:1000
    Initial_value=2*Initial_value-1;
    if Initial_value>100
        break
    end
end
params.nA   =  Initial_value; % # asset grid points (we do allow for off-grid decisions)
% put more points at lower end of asset grid. Note that optimal savings
% decisions become linear in the limit (i.e. as assets go to infinity).
% Linearly spaced grid from negative number to one
linsp_grid = linspace(linsp_min,1,params.nA);

for s=1:1000 
    Initial_value=2*Initial_value-1;
    if Initial_value>500
        break
    end
end
nA_fine=Initial_value;
% Linearly spaced grid from negative number to one
linsp_grid = linspace(linsp_min,1,params.nA);

[~, indic_bmin]=min(abs(linsp_grid));

if min(abs(linsp_grid))<10^(-8)
    linsp_grid(indic_bmin)=0;
end

% Fine Linearly spaced grid from negative number to one
linsp_grid_fine = linspace(linsp_min,1,nA_fine);
[~, indic_bmin_fine]=min(abs(linsp_grid_fine));
if min(abs(linsp_grid_fine))<10^(-8)
    linsp_grid_fine(indic_bmin_fine)=0;
end
if floor(nA_fine)~=nA_fine || linsp_grid_fine(indic_bmin_fine)~=0 || ...
        linsp_grid(indic_bmin)~=0
    error('grids don''t include borrowing constraint properly')
end
%
params.Agrid = params.bmin+(params.Amax-params.bmin)*...
    (linsp_grid.^curv);

params.Amin = min(params.Agrid);

params.Agrid = params.Agrid';
params.nA   = length(params.Agrid);
nA = params.nA;

curv_fine=curv;
params.curv_fine=curv_fine;
params.Agrid_fine = params.bmin+(params.Amax-params.bmin)*...
    (linsp_grid_fine.^curv_fine);
params.Agrid_fine = params.Agrid_fine';
nA_fine   = length(params.Agrid_fine);
params.nA_fine = nA_fine;
params.tol = 1e-12; % numerical error tolerance for consumption policy function in EGM
params.maxIter = 5*1e3; % maximum # of iterations savings problem
params.Index_b_min=find(params.Agrid==params.bmin);% Index to find in the grid bmin
params.Index_b_min_fine=find(params.Agrid_fine==params.bmin);% Index to find in the fine grid bmin

% Choose how gross liabilities and assets are spread across distribution
params.k_const_hat = 0.0;

% Select the model with flexible exchange rate and foreign denominated debt
KK = 2;

%% II.a Find good initial guess for Aiyagari equilibrium solution

if Ind_sol==0
    % Initial guesses for the parameters
    beta_vec    = (1/(1+params.r_ss))-.1:.001:(1/(1+params.r_ss))+.02; 
    N_beta      = length(beta_vec);
    resid_vec   = zeros(1,N_beta);

    for i_beta = 1:N_beta

        % Solve equilibrium for all guesses of beta
        [resid,~,~,~,~,~]=...
            calibrate_beta_open_tg([beta_vec(i_beta),1],params,NaN,NaN);
        if resid>0
            % The economy saves more than it should for low beta
            resid_vec(i_beta) = 1;
        end
        if i_beta==1 && resid_vec(1)==0
            error('Initial value of beta is too high. Try new guesses.')
        end
        if resid<0
            % Stop searching: the economy saves too little
            break
        end
    end

    % Find value of beta on the grid that is near a solution
    Ind_zz = resid_vec(2:end)-resid_vec(1:end-1);
    Ind=find(Ind_zz==-1);

    % Impose initial guess to fine solution for beta
    beta_guess = beta_vec(Ind);
    
else 
    beta_guess = 0.9832;    % Speeds up the code
end

%% II.b Solve SS by calling function that computes Aiyagari equilibrium

b_agg_ss_initial    = 0;
b_agg_ss_calibrated = params.b_agg_ss;


% Solve Steady State of the economy for beta
x0 = [beta_guess]; % Initial guess for beta derived above
options=optimset('Display','iter','MaxIter',20);
XT = fsolve(@(x) calibrate_beta_open_tg(x,params,NaN,NaN),x0,options);
[resid,dec,cpol,G0,c_fine,dec_fine,ave_mpc] = ...
    calibrate_beta_open_tg(XT,params,NaN,NaN);

params.beta = XT(1);

% This is the subjective discount factor consisten with SS targets
beta=params.beta;
beta_guess=params.beta;
save sol_ss beta_guess
if params.y_ss-sum(sum(G0.*c_fine))<0
    error('H consumption higher than output')
end

%% II.c Aggregate consumption

% Aggregate consumption of F-good and H-good
cf_agg_ss   = ((params.omega)/(1-params.omega))*...
    params.T_ss^(params.theta)*sum(c_fine(:).*G0(:));
ch_agg_ss   = sum(c_fine(:).*G0(:));

% Aggregate consumption
c_ss        = ((params.omega)^(1/params.theta)*...
    cf_agg_ss^((params.theta-1)/params.theta) +   ...
    (1-params.omega)^(1/params.theta)*...
    ch_agg_ss^((params.theta-1)/params.theta))...
    ^(params.theta/(params.theta-1));

% Fine policy function for F-good consumption
cf_fine     = ((params.omega)/(1-params.omega))*...
    params.T_ss^(params.theta)*(c_fine(:));

% SS Export of H-good
ch_star_agg = params.y_ss - sum(c_fine(:).*G0(:)) -params.k_ss*params.delta;

% Foreign demand shifter consistent with export
demshift_agg = ch_star_agg;
params.d_s =demshift_agg;

params.theta_star = ceil(-(1+(params.omega)./(1-params.omega)).^(-1)*params.b_agg_ss*(1+params.r_ss)./(params.d_s));


%% II.d Aggregate consumption
    
params.beta = beta_guess;
x0 = params.r_ss;
% Check that we have the same steady state
XT = fsolve(@(x) calibrate_R_open_tg(x,params,NaN,b_agg_ss_calibrated),x0,options);
[resid,dec,cpol,G0,c_fine,dec_fine,ave_mpc,T_calibrated,k_l_calibrated,l_calibrated,aggregates_calibrated] =...
    calibrate_R_open_tg(XT,params,NaN,b_agg_ss_calibrated);
r_calibrated=XT;

TT          = 200;          % Length of transition back to SS
rho_shock   = 0.98;           % persistence of the shock;

shock       = NaN(TT,1);    % initialize the shock vector

% Generate generic AR(1) shock
shock(1)=-params.cur_acc*params.y_ss;
for s=2:TT
    shock(s)= shock(s-1)-params.cur_acc*params.y_ss*rho_shock^(s-1);
end
b_agg_ss_final = shock(s);

% Compute the final steady state with a 2*2007 NFA
XT = fsolve(@(x) calibrate_R_open_tg(x,params,NaN,b_agg_ss_final),x0,options);
[resid,dec_final,cpol_final,G0_final,c_fine_final,dec_fine_final,ave_mpc_final,T_final,k_l_final,l_final,aggregates_final] =...
    calibrate_R_open_tg(XT,params,NaN,b_agg_ss_final);
r_final=XT;
% Compute the initial steady state with a 0 NFA
XT = fsolve(@(x) calibrate_R_open_tg(x,params,NaN,b_agg_ss_initial),x0,options);
[resid,~,~,G0_initial,~,~,ave_mpc,T_initial,k_l_initial,l_initial,aggregates_initial] =...
    calibrate_R_open_tg(XT,params,NaN,b_agg_ss_initial);
r_initial=XT;

%% 

Total_distr=G0*ones(params.ns,1);
Cum_Tot_distr=cumsum(Total_distr);
Quintile=[0.1 .3 .5 .7 .85 .95];
for j=1:length(Quintile)
   Indicator_quint(j)= min(find(Cum_Tot_distr>Quintile(j)));
   Model_nw_quintiles(j)=params.Agrid_fine(Indicator_quint(j))./params.y_ss;
end

Model_nw_quintiles
target.nw_quintiles./target.gdp

%% III.a Generate shock

% Impose AR(1) shock to external asset supply
shock_b_agg = shock;

% Shock to productivity (constant)
z_t         = ones(TT,1)*params.z_ss;


params.phi_ac   = 17; % Investment adj't cost parameter

save SteadyStates

% Technical Parameters for MIT-shock transition
% Pre-allocate endogenous variables during transition
a_prime_t   = NaN(TT,1);
c_t         = NaN(TT,1);
y_t         = NaN(TT,1);
a_t         = NaN(TT,1);
r_t_min_one = NaN(TT,1);
t_t         = NaN(TT,1);
y_t         = NaN(TT,1);

% Technical parameter - Speed of update of price guess 
lambda_par  = .005;

% III.b Transition between steady states, flexible exchange rate
% Initialize iteration
err =100;
tol = 1e-5*TT;
maxiter = 200;
params.tol_errore_t=1e-6;

for kk=KK:KK
 MIT_transition
end

save exp_flex

%% Leverage

load exp_flex
params.k_const_hat = 1/16;
k_l_exp = kl_exp;
kl_exp=fsolve(@(x) solve_transition_bis(x,z_t,new_shock_b_agg,cpol,c_fine,G0_shock,params), k_l_exp,optimset('Display','iter','UseParallel',true,'MaxFunEvals',maxfunctions));
k_l_exp = kl_exp;
 [resid,r_k_exp,r_t_exp,a_prime_exp,a_t_agg_exp,y_agg_exp,c_t_agg_exp,net_c_exp,tot_exp,c_t_exp,...
    lab_guess_exp,In_dis_exp,w_nr_guess_exp,dec_tran_exp,dec_tran_fine_exp,c_tran_fine_exp,Profits_exp,q_agg_exp]...
    =solve_transition_bis(kl_exp,z_t,new_shock_b_agg,cpol,c_fine,G0_shock,params);

k_agg_exp=[params.k_initial;kl_exp(1:TT-1)];
lab_agg_exp=kl_exp(TT:end);

save exp_flex_lev16

%% Fixed exchange rate

load exp_flex_lev16.mat

fixed=1;
params.fixed=fixed;
k_l_t_fixed=[kl_exp;tot_exp'];

[resid,r_k_fixed,r_t_fixed,a_prime_fixed,a_t_agg_fixed,y_agg_fixed,c_t_agg_fixed,net_c_fixed,tot_fixed,c_t_fixed,...
    lab_guess_fixed,In_dis_fixed,w_nr_guess_fixed,dec_tran_fixed,dec_tran_fine_fixed,c_tran_fine_fixed,Profits_fixed,q_agg_fixed]...
    =solve_transition_fixed(k_l_t_fixed,z_t,new_shock_b_agg,cpol,c_fine,G0_shock,params);

klt_fixed=fsolve(@(x) solve_transition_fixed(x,z_t,new_shock_b_agg,cpol,c_fine,G0_shock,params), k_l_t_fixed,optimset('Display','iter','UseParallel',true,'MaxFunEvals',maxfunctions));

k_l_t_fixed=klt_fixed;
[resid,r_k_fixed,r_t_fixed,a_prime_fixed,a_t_agg_fixed,y_agg_fixed,c_t_agg_fixed,net_c_fixed,tot_fixed,c_t_fixed,...
    lab_guess_fixed,In_dis_fixed,w_nr_guess_fixed,dec_tran_fixed,dec_tran_fine_fixed,c_tran_fine_fixed,Profits_fixed,q_agg_fixed]...
    =solve_transition_fixed(k_l_t_fixed,z_t,new_shock_b_agg,cpol,c_fine,G0_shock,params);

k_agg_fixed=[params.k_initial;klt_fixed(1:TT-1)];
lab_agg_fixed=klt_fixed(TT:2*TT-1);
tot_fixed=klt_fixed(2*TT:end);
savename=['exp_both_' int2str(shock_start_t) '_lev16'];
save(savename)

%% Produce figures 1-5 from the paper

Figures
