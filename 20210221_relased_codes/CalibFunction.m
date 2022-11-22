function fun = CalibFunction(x,S0,T,r,K,P)
% function fun = CalibFunction(x,S0,T,r,K,P)
% Objective function for market calibration of Variance Gamma Model

theta = x(1);
nu = x(2);
sigma = x(3);

% Risk-neutral Parameter
omega = (1/nu).*...
    log(1-theta*nu - sigma*sigma*nu/2);

% Get distinct Maturities
Tdistinct = unique(T);
Nmaturities = length(Tdistinct);

fun = nan(numel(K),1);

for i = 1:Nmaturities
    % Get maturities
    idx = T == Tdistinct(i);
    
    
    % perform pricing via FFT
    [CallPrices,LogStrikes] = FFTPricing(Tdistinct(i),r,...
        @(w)phi_vg(w,S0,r,omega,Tdistinct(i),theta,nu,sigma));
    
    % Get call Prices corrisponding to given prices
    K_sel = K(idx);
    
    % Prezzo delle Call
    Call_Prices = interp1(LogStrikes,CallPrices,log(K_sel));
    
    % Funzione da minimizzare
    % fun = (Call_Prices - P(idx)).^2; (this is wrong)
    fun(idx) = (Call_Prices - P(idx)).^2;
    
end




end