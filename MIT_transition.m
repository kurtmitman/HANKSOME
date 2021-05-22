%% MIT_transition
% computes the transition path from the initial to the final steady state
% as well as the unexpected contraction in foreign credit

kk
 
if kk==1
    fixed=false;
    foreign=false;
    string_name{kk}='Flexible, domestic';
elseif kk==2
    fixed=false; 
    foreign=true;
    string_name{kk}='Flexible, foreign';
elseif kk==3
    fixed=true; 
    foreign=true;
    string_name{kk}='Fixed, foreign';
end
    
lambda       = lambda_par*linspace(1,0.01,TT).^.9';
iter_mkt_clr = 1;

r_k_guess = ones(TT,1);
t_guess = ones(TT,1);    
k_flex_guess = ones(TT,1);
l_flex_guess = ones(TT,1);
for t=1:TT
    r_k_guess(t)    = r_final+(r_initial-r_final)*0.8^t;
    t_guess(t)      = T_final+(1.1-T_final)*0.8^(t-1);
    k_flex_guess(t) = k_l_final*l_final+(k_l_initial*l_initial-k_l_final*l_final)*0.9^(t);
    l_flex_guess(t) = l_final+(l_initial-l_final)*0.9^(t-1);
end

% Preallocate generic-period policy functions during transition
c_new           = NaN(nA,nx);
c_new_fine      = zeros(nA_fine,ns);
dec_temp_fine   = zeros(nA_fine,ns);

% Preallocate policy functions for whole transition
c_tran          = NaN(nA,nx,TT);
c_tran_fine     = NaN(nA_fine,nx,TT);
dec_tran        = NaN(nA,nx,TT);

% Impose SS policy functions in last period of transition
c_tran(:,:,TT)      = cpol_final;
c_tran_fine(:,:,TT) = c_fine_final;

% Repmat savings vector to state-space-sized grid
sav_grid    = repmat(params.Agrid,[1 ns]);

% Repmat exog.state vector to state-space-sized grid
exLw        = repmat(params.ex',[nA 1]);

% Grid and fine grid for assets
Agrid       = params.Agrid;
Agrid_fine  = params.Agrid_fine;

% Grids of gross assets and liabilities
Bgrid_fine =NaN(size(Agrid_fine));
Kgrid_fine =NaN(size(Agrid_fine));

Bgrid_fine(Agrid_fine<0)=Agrid_fine(Agrid_fine<0);
Kgrid_fine(Agrid_fine<0)=0;

% Total Wealth of the positive wealth holders
net_a_rich = (sum(G0(Agrid_fine>0,:)')*Agrid_fine(Agrid_fine>0));

% Leverage ratio of the rich
lev_rat_rich = params.k_ss/net_a_rich;

% Gross assets of the rich
Kgrid_fine(Agrid_fine>=0)=lev_rat_rich*Agrid_fine(Agrid_fine>=0);

% Gross liabilities of the rich
Bgrid_fine(Agrid_fine>=0)=-(Kgrid_fine(Agrid_fine>=0)-Agrid_fine(Agrid_fine>=0));

%% Transition from one SS to another
iter_mkt_clr=1;
close all
 
% Preallocate grid to compute new distribution
A_grid_fine_new=NaN(size(Agrid_fine,1),size(Agrid_fine,2));
err=100;
params.fixed=fixed;
params.TT=TT;
params.foreign=foreign;
params.linsp_min = linsp_min;
params.k_initial = k_l_initial*l_initial;
params.k_final = k_l_final*l_final;
params.i_ss = params.delta*params.k_initial;

load transition_start.mat
params.BG_shock=0;
params.q_initial=1;
params.tot_initial=T_initial;
k_l = kl_path;

if params.fixed==false
    [resid,r_k,r_t,a_prime,a_t_agg,y_agg,c_t_agg,net_c,tot,c_t,...
    lab_guess,In_dis,w_nr_guess,dec_tran,dec_tran_fine,c_tran_fine,Profits,q_agg,dist_tran]...
    =solve_transition_bis(k_l,z_t,shock_b_agg,cpol_final,c_fine_final,G0_initial,params);
if max(abs(resid))>10^(-6)
    kl_path=fsolve(@(x) solve_transition_bis(x,z_t,shock_b_agg,cpol_final,c_fine_final,G0_initial,params), k_l,optimset('Display','iter','UseParallel',true));
    [resid,r_k,r_t,a_prime,a_t_agg,y_agg,c_t_agg,net_c,tot,c_t,...
    lab_guess,In_dis,w_nr_guess,dec_tran,dec_tran_fine,c_tran_fine,Profits,q_agg,dist_tran]...
    =solve_transition_bis(kl_path,z_t,shock_b_agg,cpol_final,c_fine_final,G0_initial,params);
else
   kl_path=k_l;
end

else
    [resid,r_k,r_t,a_prime,a_t_agg,y_agg,c_t_agg,net_c,tot,c_t,...
    lab_agg,In_dis,w_r,dec_tran,dec_tran_fine,c_tran_fine]...
    =solve_transition_fixed(k_l,z_t,shock_b_agg,cpol_final,c_fine_final,G0_initial,params);
if max(abs(resid))>10^(-6)
    kl_path=fsolve(@(x) solve_transition_fixed(x,z_t,shock_b_agg,cpol_final,c_fine_final,G0_initial,params), k_l,optimset('Display','iter','UseParallel',true));
    [resid,r_k,r_t,a_prime,a_t_agg,y_agg,c_t_agg,net_c,tot,c_t,...
    lab_guess,In_dis,w_nr_guess,dec_tran,dec_tran_fine,c_tran_fine]...
    =solve_transition_fixed(kl_path,z_t,shock_b_agg,cpol_final,c_fine_final,G0_initial,params);
else
    kl_path=k_l;
end
end    
save transition_start_new kl_path

%% Unexpected credit contraction
k_agg=[aggregates_initial.k;kl_path(1:TT-1)];
lab_agg=kl_path(TT:end);
shock_start_t = 41;
temp = linspace(1,0,TT)';
new_shock_b_agg = aggregates_calibrated.BGss+(shock_b_agg(shock_start_t-1)-aggregates_calibrated.BGss)*temp.^4;

k_exp_guess = ones(TT,1);
l_exp_guess = ones(TT,1);
for t=1:TT
    k_exp_guess(t) = aggregates_calibrated.k+(k_agg(shock_start_t)-aggregates_calibrated.k)*0.8^(t);
    l_exp_guess(t) = aggregates_calibrated.l;
end

k_l_exp = [k_exp_guess(1:TT-1);l_exp_guess];
G0_shock=squeeze(dist_tran(:,:,shock_start_t));
params.k_initial = k_agg(shock_start_t);
params.k_final = aggregates_calibrated.k;
params.i_ss = params.delta*params.k_final;
params.BG_shock=shock_b_agg(shock_start_t-1);
params.q_initial=q_agg(shock_start_t-1);
params.tot_initial=tot(shock_start_t-1);

maxfunctions=5000;
kl_exp=fsolve(@(x) solve_transition_bis(x,z_t,new_shock_b_agg,cpol,c_fine,G0_shock,params), k_l_exp,optimset('Display','iter','UseParallel',true,'MaxFunEvals',maxfunctions));
 [resid,r_k_exp,r_t_exp,a_prime_exp,a_t_agg_exp,y_agg_exp,c_t_agg_exp,net_c_exp,tot_exp,c_t_exp,...
    lab_guess_exp,In_dis_exp,w_nr_guess_exp,dec_tran_exp,dec_tran_fine_exp,c_tran_fine_exp,Profits_exp,q_agg_exp]...
    =solve_transition_bis(kl_exp,z_t,new_shock_b_agg,cpol,c_fine,G0_shock,params);

k_agg_exp=[params.k_initial;kl_exp(1:TT-1)];
lab_agg_exp=kl_exp(TT:end);

