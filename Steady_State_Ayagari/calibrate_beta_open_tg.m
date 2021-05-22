function [resid,dec,cpol,G0,c_fine,dec_fine,ave_mpc] =...
    calibrate_beta_open_tg(XT,params,c0,G0)
    % Two good open economy with heterogeneous agents


    Agrid = params.Agrid;   % Asset grid
    ns = params.ns;         % Number of exogenous state points
    nA = params.nA;         % Number of endogenous state points
    ex = params.ex;         % Exogenous state values' vector
    exinv = params.exinv;   % Exogenous state ergodic distribution

    nA_fine = params.nA_fine;   % Number of end. state fine-grid points
    Agrid_fine = params.Agrid_fine;     % End. state fine-grid

    exL = ex; 

    beta = XT(1);

    params.beta = beta;

    BGss = params.b_agg_ss; % SS external assets
    Kss  = params.k_ss;  % SS capital stock

    r    = params.r_ss;  % SS interest rate
    prof = params.pi_ss;  % SS firms' profits

    % Initial guess for consumption policy function
    if(isnan(c0)==1)
        disp('No initial guess for c0')
        c0 = (1+params.T_ss^(params.theta-1)*...
            (params.omega/(1-params.omega)))^(-1)...
            *(params.w_ss *repmat(exL',[nA 1]) +...
            r*repmat(Agrid,[1 ns]) + prof);
    end
    
    % stochastic consumption savings problem
    [dec,cpol,dec_fine,c_fine] = solve_EGM_EL_open_tg(r,exL,beta,c0,params); 

    ComputeDistHist_open   

    A = sum(a);
    CH = sum(c);
    ave_mpc = sum(mpc);
    beta
    
    resid(1)   = BGss + Kss - A;

resid
end