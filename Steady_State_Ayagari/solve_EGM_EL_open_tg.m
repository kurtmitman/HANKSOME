function [dec,cpol,dec_fine,c_fine] = solve_EGM_EL_open_tg(r,exL,beta,c0,params)
% solve_EGM_EL_open_tg Function that solves stochastic consumption savings problem
% using the endogenous grid point method. Two-good open economy.
%
% INPUT:        w     wage rate (NO)
%               r     interest rate
%               exL   labor endoemwnt
%               c0    initial guess consumption function
%               param parameter structure
%
% OUTPUT:       dec   optimal savings decisions (nA*ns)
%               cpol  optimal T-consumption decisions (nA*ns)
%               dec_fine same as above, but for fine grid (nA_fine*ns)
%               c_fine same as above, but for fine grid (nA_fine*ns)


nA = params.nA;
nA_fine = params.nA_fine;
ns = params.ns; 
Agrid = params.Agrid;
Agrid_fine = params.Agrid_fine;
tol = params.tol;
maxIter = params.maxIter;

sav_grid = repmat(Agrid,[1 ns]);

exLw = repmat(exL',[nA 1]);

c_constrained =max((1+params.T_ss^(params.theta-1)*(params.omega/(1-params.omega)))^(-1)*( params.w_ss*exLw*params.l_bar ...
    +(1+params.r_ss)*repmat(Agrid,[1 ns]) + params.pi_ss - params.bmin),10^(-5));

for m=1:maxIter

    Emup1 = c0.^(-params.gamma);
    Emup = beta*(1+params.r_ss)*Emup1*params.piex';
    c_s = Emup.^(-1/params.gamma);
    a_today = (c_s*(1+params.T_ss^(params.theta-1)*(params.omega/(1-params.omega)))+sav_grid - params.w_ss*exLw*params.l_bar  - params.pi_ss)/(1+r);
    
    
    for j=1:ns
        
        c_new(:,j) = (Agrid>a_today(params.Index_b_min,j)).*interp1(a_today(:,j),c_s(:,j),Agrid,'pchip')+(Agrid<=a_today(params.Index_b_min,j)).*c_constrained(:,j);
        c_new(:,j)=  max(c_new(:,j),10^-5);
    end
    
    maxDiff = max(max(abs(c_new-c0)));
    if(maxDiff < 1e-4)
        updateval = 0.5;
    else
        updateval = 0.9;
    end
    
    if (maxDiff < tol)
        fprintf('consumption policy function converged after %d iterations \n',m);
        break
    end
    if (m==maxIter)
        fprintf('consumption policy function did NOT converge after %d iterations \n',maxIter);
    end
    c0 = c_new*updateval+(1-updateval)*c0;

end

dec  = (1+r)*repmat(Agrid, [1 ns])+exLw*params.w_ss*params.l_bar...
    -c0*(1+params.T_ss^(params.theta-1)*(params.omega/(1-params.omega))) + params.pi_ss;
cpol = c0;
c_fine = zeros(nA_fine,ns);
dec_fine = zeros(nA_fine,ns);
for j=1:ns
    c_fine(:,j) = (Agrid_fine>a_today(params.Index_b_min,j)).*interp1(a_today(:,j),c_s(:,j),Agrid_fine,'pchip')+(Agrid_fine<=a_today(params.Index_b_min,j)).*interp1(Agrid,c_constrained(:,j),Agrid_fine,'pchip');
    dec_fine(:,j) = (Agrid_fine>a_today(params.Index_b_min,j)).*((1+r)*Agrid_fine+exL(j)...
        *params.w_ss*params.l_bar-c_fine(:,j)*(1+params.T_ss^(params.theta-1)*(params.omega/(1-params.omega)))...
        + params.pi_ss)+(params.Agrid_fine<=a_today(params.Index_b_min,j)).*params.bmin;
end


end