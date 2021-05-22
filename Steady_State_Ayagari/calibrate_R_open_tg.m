function [resid,dec,cpol,G0,c_fine,dec_fine,ave_mpc,T,k_l,l,aggregates] =...
    calibrate_R_open_tg(XT,params,c0,BGss)
    % Two good open economy with heterogeneous agents

    Agrid = params.Agrid;   % Asset grid
    ns = params.ns;         % Number of exogenous state points
    nA = params.nA;         % Number of endogenous state points
    ex = params.ex;         % Exogenous state values' vector
    exinv = params.exinv;   % Exogenous state ergodic distribution

    nA_fine = params.nA_fine;   % Number of end. state fine-grid points
    Agrid_fine = params.Agrid_fine;     % End. state fine-grid

    exL = ex; 

r    = XT;  % SS interest rate
r_K  =r+params.delta; % SS return on capital, net of depreciation
k_l=(r_K./params.alpha.*(params.epsilon_f)/(params.epsilon_f-1)).^(1/(params.alpha-1));
w = k_l.^params.alpha.*(1-params.alpha).*(params.epsilon_f-1)/(params.epsilon_f);
  
rhs_ch=@(x)(params.epsilon_w-1)/params.epsilon_w*(1-params.omega)*w ./params.chi_dis.* x^(-params.varphi);
lhs_ch=@(x)(1+params.omega/(1-params.omega)).^(-1).*(r*(BGss+k_l.*x)...
    +params.z_ss*x*k_l^params.alpha.*((params.epsilon_f -1)./params.epsilon_f.*(1-params.alpha)+1/params.epsilon_f));
find_l=@(x)[rhs_ch(x)-lhs_ch(x)];
x=fsolve(find_l,1);
l=x;
c_h=rhs_ch(x);
k=k_l*l;
y=params.z_ss*k^(params.alpha)*l^(1-params.alpha);
Profits=1/params.epsilon_f*y;
T=((y - c_h -k*params.delta)/params.d_s).^(-1/params.theta_star);  

    % Initial guess for consumption policy function
    if(isnan(c0)==1)
        disp('No initial guess for c0')
        c0 = (1+T^(params.theta-1)*...
            (params.omega/(1-params.omega)))^(-1)...
            *(w*l *repmat(exL',[nA 1]) +...
            r*repmat(Agrid,[1 ns]) + Profits);
    end
    
    % stochastic consumption savings problem
    [dec,cpol,dec_fine,c_fine] = solve_R(r,exL,l,w,T,Profits,c0,params);
    ComputeDistHist_open   


    A = sum(a);
    CH = sum(c);
    ave_mpc = sum(mpc);
    
    
    resid(1)   = BGss + k - A;
    aggregates.Profits=Profits;
    aggregates.y=y;
    aggregates.k=k;
    aggregates.c_h=c_h;
    aggregates.l=l;
    aggregates.w=w;
    aggregates.r_k=r_K;
    aggregates.BGss=BGss;

resid
end