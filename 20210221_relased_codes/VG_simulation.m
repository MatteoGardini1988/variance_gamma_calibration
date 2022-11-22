function X = VG_simulation(Nsim,nDates,dt,params)
theta = params(1);
nu = params(2);
sigma = params(3);
g = gamrnd(dt/nu,nu,Nsim,nDates);

% S = S0*cumprod([ones(Nsim,1) exp(mu*dt + theta*g + sigma*sqrt(g).*randn(Nsim,T))],2);
% X = cumsum([zeros(Nsim,1) theta*g + sigma*sqrt(g).*randn(Nsim,T)],2);

X = zeros(Nsim,nDates);

for j = 2:nDates
    X(:,j) = X(:,j-1) + theta.*g(:,j) + sigma.*sqrt(g(:,j)).*randn(Nsim,1); 
end


end