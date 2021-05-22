function [c_t,a_prime_t,dec_tran,dec_tran_fine,c_tran_fine,In_dis,dist_tran]= backsolve_egm(r_t,t_guess,y_agg,cpol,c_fine,prof_firm,Profit_i,G0,params);

nA          = params.nA;
nA_fine     = params.nA_fine;
nx          = params.nx;
ns          = params.ns;
TT = params.TT;
beta        = params.beta;
foreign     = params.foreign;
linsp_min   = params.linsp_min;

% Preallocate generic-period policy functions during transition
c_new           = NaN(nA,nx);
c_new_fine      = zeros(nA_fine,ns);
dec_temp_fine   = zeros(nA_fine,ns);

% Preallocate policy functions for whole transition
c_tran          = NaN(nA,nx,TT);
c_tran_fine     = NaN(nA_fine,nx,TT);
dec_tran        = NaN(nA,nx,TT);

% Impose SS policy functions in last period of transition
c_tran(:,:,TT)      = cpol;
c_tran_fine(:,:,TT) = c_fine;

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

params.k_const = params.k_const_hat*params.k_initial;

% Total Wealth of the positive wealth holders
net_a_rich = (sum(G0(Agrid_fine>0,:)')*Agrid_fine(Agrid_fine>0));

% Leverage ratio of the rich, for variable capital
lev_rat_rich = (1-params.k_const_hat)*params.k_initial*params.q_initial/net_a_rich;
if(params.BG_shock==0)
    lev_rat_rich=1;
end
    
% Gross assets 
Kgrid_fine = max(lev_rat_rich*Agrid_fine + params.k_const,0);

% Gross liabilities 
Bgrid_fine = -(Kgrid_fine - Agrid_fine);


for t=TT-1:-1:1
    
    % Interest rate between t and t+1
    r           =   r_t(t+1);
    % Interest rate between t-1 and t
    r_yesterday =   r_t(t);
    tot     =t_guess(t);
    y_t_agg =y_agg(t); 
    prof_firms_t=prof_firm(t)+Profit_i(t);

    % Consumption of agent against borrowing constraint
    c_constrained =max( (1+tot.^(params.theta-1)*(params.omega/(1-params.omega))).^(-1)*(...
        (1+r_yesterday)*repmat(Agrid,[1 ns]) +(params.epsilon_f-1)./params.epsilon_f.* exLw*y_t_agg*(1-params.alpha)/params.ex_mean...
        +prof_firms_t - params.bmin),10^(-5));
    
    Emup1 = c_tran(:,:,t+1).^(-params.gamma);   % Future MU of consumption
    Emup = beta*(1+r)*Emup1*params.piex';       % Expected MU of fut. cons.
    c_s = Emup.^(-1/params.gamma);              % Current cons given future assets
    
    % Current assets given current cons and future assets (budget const.)
    a_today = (c_s*(1+tot.^(params.theta-1)*(params.omega/(1-params.omega)))+sav_grid -(params.epsilon_f-1)./params.epsilon_f*exLw*y_t_agg*(1-params.alpha)/params.ex_mean...
        -prof_firms_t)/(1+r_yesterday);
    
    % Update consumtpion policy function via interpolation
    for j=1:ns
        
        c_new(:,j) = (Agrid>a_today(params.Index_b_min,j)).*interp1(a_today(:,j),c_s(:,j),Agrid,'pchip')...
            + (params.Agrid<=a_today(params.Index_b_min,j)).*c_constrained(:,j);
        c_new(:,j)=  max(c_new(:,j),10^-5);
    end
    
    % Allocate period-t consumption p.f.
    c_tran(:,:,t) = c_new;
    
    % Allocate period-t p.f. for next-period assets
    dec_tran(:,:,t) =  (1+r_yesterday)*repmat(Agrid, [1 ns])+...
        (params.epsilon_f-1)./params.epsilon_f*exLw*y_t_agg*(1-params.alpha)/params.ex_mean...
        +prof_firms_t -c_new.*(1+tot.^(params.theta-1)*(params.omega/(1-params.omega)));
    
    % Update fine p.f.'s for consumption and saving
    for j=1:ns
        
        c_new_fine(:,j)     = (Agrid_fine>a_today(params.Index_b_min,j)).*interp1(a_today(:,j),c_s(:,j),Agrid_fine,'pchip')...
            +(params.Agrid_fine<=a_today(params.Index_b_min,j)).*interp1(Agrid,c_constrained(:,j),Agrid_fine,'pchip');
        
        dec_temp_fine(:,j)  = (Agrid_fine>a_today(params.Index_b_min,j)).*...
            ((1+r_yesterday)*Agrid_fine +(params.epsilon_f-1)./params.epsilon_f* params.ex(j)*y_t_agg*(1-params.alpha)/params.ex_mean...
            +prof_firms_t ...
            - c_new_fine(:,j)*(1+tot.^(params.theta-1)*(params.omega/(1-params.omega)) ))...
            +(params.Agrid_fine<=a_today(params.Index_b_min,j)).*params.bmin;
        
    end
   if  min(dec_temp_fine(:))<params.bmin
       keyboard
       error('BEEEP')
   end
    % Allocate period-t fine consumption and assets p.f.'s
    c_tran_fine(:,:,t)      = c_new_fine;
    dec_tran_fine(:,:,t)    = dec_temp_fine;
  
end

%%
% Compute evolution of asset distribution
% Preallocate distribution during transition
dist_tran(:,:,1)    = G0;
if foreign

    % Debt re-evaluation
    Bgrid_fine_new = Bgrid_fine/(t_guess(1)/params.tot_initial);

    % New net assets after revaluation
    A_grid_fine_new = Bgrid_fine_new + Kgrid_fine;
   
    if params.bmin/t_guess(1)<Agrid_fine(1)
        error('There is positive mass in points outside Agrid_fine')
    end

    G0_new=0*G0;

    for ai = 1:nA_fine
        [vals,inds]=basefun(Agrid_fine,nA_fine,A_grid_fine_new(ai));
        for j=1:ns
          G0_new(inds(1),j) = G0_new(inds(1),j)+G0(ai,j)*vals(1);
          G0_new(inds(2),j) = G0_new(inds(2),j)+G0(ai,j)*vals(2);
        end
    end
    
    if t_guess(1)>1
             G0_new( find(isnan(G0_new)==1))=0;
    end
    if sum(G0_new(:))<1-1e-6
        error('The new distribution does not sum to one')
    end
    dist_tran(:,:,1)    = G0_new;
end

%% We need to change the transition

dist_tran(:,:,2:TT) = 0;

% Initial period of transition. Compute aggregate saving and consumption
t=1;
a_prime = sum(dist_tran(:,:,t).*dec_tran_fine(:,:,t),1);
c       = sum(dist_tran(:,:,t).*c_tran_fine(:,:,t),1);

% Store initial-period aggregate saving and consumption
a_prime_t(t) = sum(a_prime);
c_t(t) = sum(c);

%% UNDERSTAND VALA!
for t=2:TT-1   
    
   for j=1:ns
        apl = floor( (nA_fine-1)*((( dec_tran_fine(:,j,t-1)-params.bmin)/...
            (params.Amax-params.bmin)).^(1/params.curv_fine) - linsp_min )/...
            (1-linsp_min))+1;
        apl=apl.*(apl<nA_fine)+(nA_fine-1)*(apl>=nA_fine);
        vala = (params.Agrid_fine(apl+1)- dec_tran_fine(:,j,t-1))./...
            (params.Agrid_fine(apl+1)- params.Agrid_fine(apl));
        for ai = 1:nA_fine     
            for jp=1:ns
              dist_tran(apl(ai),jp,t) = dist_tran(apl(ai),jp,t)+...
                  dist_tran(ai,j,t-1)*params.piex(j,jp)*vala(ai);
              dist_tran(apl(ai)+1,jp,t) = dist_tran(apl(ai)+1,jp,t)+...
                  dist_tran(ai,j,t-1)*params.piex(j,jp)*(1-vala(ai));
                
            end
        end
   end

     % Store aggregate saving and consumption   
     a_prime    = sum(dist_tran(:,:,t).*dec_tran_fine(:,:,t),1);
     c          = sum(dist_tran(:,:,t).*c_tran_fine(:,:,t),1);
     a_prime_t(t)   = sum(a_prime);
     c_t(t)         = sum(c);

In_dis=dist_tran(:,:,1);
end
