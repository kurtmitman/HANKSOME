function [resid,r_k,r_t,a_prime,a_t_agg,y_agg,c_t_agg,net_c,tot,c_t,...
    lab_agg,In_dis,w_r,dec_tran,dec_tran_fine,c_tran_fine,Profits,q_agg,dist_tran]...
    =solve_transition_fixed(XX,z_t,shock_b_agg,cpol,c_fine,G0,params)

foreign=params.foreign;
k_agg = [params.k_initial;XX(1:params.TT-1)];
lab_agg=XX(params.TT:2*params.TT-1);
t_guess=XX(2*params.TT:3*params.TT-1);

TT = params.TT;
resid=zeros(2*TT-1,1);
delta = params.delta;
i_ss  = params.i_ss;
a_prime=zeros(TT,1);
a_t_agg=zeros(TT,1);
i_agg = k_agg(2:TT)- (1-delta) *k_agg(1:TT-1);
q_agg =[ 1+ params.phi_ac*(i_agg(1:TT-1)-params.delta*k_agg(1:TT-1)).*(k_agg(1:TT-1).^(-1)) ; 1];

y_agg =z_t.*(k_agg).^params.alpha.*lab_agg.^(1-params.alpha);
a_prime(1:TT-1) =  shock_b_agg(1:TT-1) + q_agg(1:TT-1).*k_agg(2:TT);
a_prime(TT) = q_agg(TT)*params.k_final + shock_b_agg(TT);

a_t_agg(1) =params.q_initial*k_agg(1)+params.BG_shock;
a_t_agg(2:TT) = a_prime(1:TT-1);

Profits_i =[(q_agg(1:TT-1)-1).*i_agg(1:TT-1) - params.phi_ac/2.*k_agg(1:TT-1).*((i_agg(1:TT-1)-delta*k_agg(1:TT-1))./k_agg(1:TT-1)).^2;
    0];


c_t_agg = zeros(TT,1);
net_c = zeros(TT,1);

        % 1 start with a guess of t
        % 2 Find the inflation implied
 
        t_inf=[t_guess(1)/params.tot_initial; t_guess(2:end)./t_guess(1:end-1)];
        % Find the marginal cost from NKPC
        mc(TT,1)=(params.epsilon_f-1)./params.epsilon_f...
            + params.theta_nr./params.epsilon_f.*(t_inf(TT)-1).*t_inf(TT);
        % Find the real wage
        w_guess(TT,1)=  mc(TT).*(params.alpha).^(params.alpha).*(1-params.alpha).^(1-params.alpha)...
            .*(params.alpha./(1-params.alpha)).^(-params.alpha).*(k_agg(TT)./lab_agg(TT)).^(params.alpha);
        % Find the MPK capital
        MPK_guess(TT,1)= params.alpha.*w_guess(TT).*lab_agg(TT)./((1-params.alpha).*k_agg(TT));
        % Find the real interest rate
        r_t(TT,1)=( MPK_guess(TT) + (1- delta)*q_agg(TT))./q_agg(TT-1)-1;
        
        for t=TT-1:-1:1
     mc(t,1)=(params.epsilon_f-1)./params.epsilon_f+ params.theta_nr./params.epsilon_f.*(t_inf(t)-1).*t_inf(t)...
         -1/(1+r_t(t+1,1)).*params.theta_nr./params.epsilon_f.*(t_inf(t+1)-1).*t_inf(t+1).*y_agg(t+1)./y_agg(t);
        w_guess(t,1)=  mc(t).*(params.alpha).^(params.alpha).*(1-params.alpha).^(1-params.alpha)...
            .*(params.alpha./(1-params.alpha)).^(-params.alpha).*(k_agg(t)./lab_agg(t)).^(params.alpha);
        MPK_guess(t,1)= params.alpha.*w_guess(t).*lab_agg(t)./((1-params.alpha).*k_agg(t));
        if t>1
            r_t(t)=( MPK_guess(t) + (1- delta)*q_agg(t))./q_agg(t-1)-1;
        else
             r_t(t)=( MPK_guess(t) + (1- delta)*q_agg(t))/params.q_initial-1;
      
        end
            
        end    
       Profits_guess=y_agg- MPK_guess.*k_agg-w_guess.*lab_agg-params.theta_nr./2.*(t_inf-1).^2;
  
            net_c(1:TT-1)=y_agg(1:TT-1)-i_agg(1:TT-1)-t_guess(1:TT-1).^(-params.theta_star)*params.d_s ...
                - params.phi_ac/2.*k_agg(1:TT-1).*((i_agg(1:TT-1)-delta*k_agg(t))/k_agg(t)).^2 - -params.theta_nr./2.*(t_inf(1:TT-1)-1).^2;
            net_c(TT)=y_agg(TT)-i_agg(TT-1)-t_guess(TT).^(-params.theta_star)*params.d_s ...
                - params.phi_ac/2*k_agg(TT)*((i_agg(TT-1)-delta*k_agg(TT))/k_agg(TT)).^2 - -params.theta_nr./2.*(t_inf(TT)-1).^2;

r_k= MPK_guess;   
w_r=w_guess;
Profits_firms=Profits_guess;


tot = t_guess;
[c_t,a_prime_t,dec_tran,dec_tran_fine,c_tran_fine,In_dis,dist_tran]...
    = backsolve_egm(r_t,tot,y_agg,cpol,c_fine,Profits_firms,Profits_i,G0,params);

c_t=[c_t';c_t(end)];
k_prime_t = (a_prime_t' - shock_b_agg(1:TT-1))./q_agg(1:TT-1);
lab_prime=(1-params.omega).*(params.epsilon_w-1)./...
            (params.epsilon_w).*w_r./(c_t.* params.chi_dis);

Profits=Profits_firms+Profits_i;
resid(1:params.TT-1)=k_prime_t-k_agg(2:params.TT);
resid(params.TT:2*TT-1)=lab_prime-lab_agg;
resid(2*params.TT:3*TT-1)=net_c-c_t;


