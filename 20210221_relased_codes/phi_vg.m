function v = phi_vg(u,S0,r,omega,T,theta,nu,sigma)


    % tmp = 1 - 1i * theta * nu * u + 0.5 * sigma * sigma * u .* u * nu;
    % y = 1i * u * (lnS + (r + omega - d) * T ) - T*log(tmp)/nu;

% GARDINI
%     v = (S0.*exp((r + omega).*T)).*...
%         (1 - 1i * theta *nu *u + 0.5.*sigma^2.*u.^2.*nu).^(-T/nu);

    argomento = 1 - 1i*theta * nu * u + 0.5 * sigma * sigma .* u .* u *nu; 
    v = exp((log(S0) + (r + omega) * T).*1i.*u) .* argomento.^(-T/nu);



end