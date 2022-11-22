function [CallPrices,LogStrikes] = FFTPricing(T,r,phi)
% [CallPrices,Strikes] = FFTPricing(S0,K,T,phi)
% Prezza delle Call Europee tramite metodo della FFT di Carr-Madan 1999
% T: maturity
% r: risk free 
% phi: funzione caratteristica del modello di pricing scelto.
% Funzione caratteristica di Black e Scholes
% phi=@(x) exp(1i*(log(S0) + (r-0.5*sigma^2)*T).*x -0.5*sigma^2*T.*x.^2);
%
% MG: 02marzo2017

% Intero (Serve per definire il numero di nodi di quadratura)
L=12;

% N potenza del due
N=2^L;

% spazio sulla griglia di integrazione
dv=0.25;
% dk Spazio sulla griglia del log(K)
dk=2*pi/(N*dv);

% Alpha
alpha=0.75;

% Vettore v (integranda)
jvec=1:N;
v= (jvec-1)*dv;

% vettore k (log strike)
uvec=1:N;

% Estremo inferiore di ku (va scelto di modo che la griglia sui k si
% concentri attorno alle opzioni at the money perchè sono quelle che più
% interessa prezzare).
b=(N * dk) / 2; 
ku=-b + dk*(uvec-1);

% Da trasformare
psi= exp(-r*T) .*phi(v-(alpha+1)*1i)./(alpha^2+alpha - v.^2 +1i*(2*alpha+1).*v);

% Quadratura di Simpson
SimpsonW=zeros(1,N) ;
SimpsonW(1)=1/3;
for i =2:N
    
    if(mod( i ,2)==0)
        SimpsonW( i )=4/3;
    else
        SimpsonW( i )=2/3;
    end
end
tmp=dv.*psi .* exp(1i * v *b ).*SimpsonW;  

% tmp=dv*(psi .* exp(1i * v *b ) - ...
%     0.5*(...
%     exp(-1i*v(1))*phi(v(1)) + ...
%     exp(-1i*v(end))*phi(v(end))));  



% Prezzo della Call
CallPrices=exp(-alpha.*ku)./pi .* real(fft(tmp));
LogStrikes=ku;


end