function[op_out, Toep, alpha_LDR] = ...
    getOPERATOR_DDA(r, ko, refInd, kvec, Eo, nearby_quad)

%% Self-interaction Classius-Mossotti stuff
[L,M,N,~]=size(r);
dx=r(2,1,1,1)-r(1,1,1,1);
b1 = -1.8915316;
b2 = 0.1648469;
b3 = -1.7700004;
msqr = refInd.^2;
dcube = dx^3;
d=dx;

% if nargin > 3  % we have polarization info   
  a_hat = kvec/norm(kvec);
  e_hat = Eo/norm(Eo);
  S = 0;
  for j = 1:3
    S = S + (a_hat(j)*e_hat(j))^2;
  end
% else           % use randomly-oriented value; also for non-plane wave
%   S = .2;
% end    

alpha_CM = 3/(4*pi)*(msqr - 1)./(msqr + 2); % Clausius-Mossotti
alpha_LDR = alpha_CM./(1 + (alpha_CM).*((b1+msqr*b2+msqr*b3*S)*(ko*d)^2-2/3*1i*ko^3*dcube));


I = eye(3);

Toep = zeros(L,M,N,3);
R0 = squeeze(r(1,1,1,:));


nQuad=10;
[wG,xG] = gauss_1d(nQuad);
[XG,YG,ZG] = meshgrid(xG);
[XW,YW,ZW] = meshgrid(wG*0.5);

for i=1:L
    for j=1:M
        for k=1:N
            R1 = squeeze(r(i,j,k,:));
            rk_to_rj = R1-R0;
            rjk = norm(rk_to_rj);
            
            if strcmp(nearby_quad, 'on') == 1
                if rjk<5*dx && rjk>1e-15
                    
                    x_grid = R1(1) + dx/2 * XG;
                    y_grid = R1(2) + dx/2 * YG;
                    z_grid = R1(3) + dx/2 * ZG;
                    
                    temp=zeros(3,3);
                    for iQ = 1:nQuad
                        for jQ = 1:nQuad
                            for kQ = 1:nQuad
                                RQ = [x_grid(iQ,jQ,kQ);y_grid(iQ,jQ,kQ);...
                                    z_grid(iQ,jQ,kQ)];
                                
                                rk_to_rj = RQ-R0;
                                
                                rjk = norm(rk_to_rj);
                                rjk_hat = (rk_to_rj)/rjk;
                                rjkrjk = rjk_hat*rjk_hat';
                                
                                Ajk = exp(1i*ko*rjk)/rjk*(ko^2*(I - rjkrjk) + (1i*ko*rjk-1)/rjk^2*(I - 3*rjkrjk)); %Draine & Flatau
                                temp = temp + Ajk.*XW(iQ,jQ,kQ).*...
                                    YW(iQ,jQ,kQ).*ZW(iQ,jQ,kQ);
                                
                            end
                        end
                    end
                    Toep(i,j,k,1) = temp(1,1);
                    Toep(i,j,k,2) = temp(1,2);
                    Toep(i,j,k,3) = temp(1,3);
                    Toep(i,j,k,4) = temp(2,2);
                    Toep(i,j,k,5) = temp(2,3);
                    Toep(i,j,k,6) = temp(3,3);
                else
                    rk_to_rj = R1-R0;
                    
                    rjk = norm(rk_to_rj);
                    rjk_hat = (rk_to_rj)/rjk;
                    rjkrjk = rjk_hat*rjk_hat';
                    
                    if abs(rjk)>1e-15
                        Ajk = exp(1i*ko*rjk)/rjk*(ko^2*(I - rjkrjk) + (1i*ko*rjk-1)/rjk^2*(I - 3*rjkrjk)); %Draine & Flatau
                        Toep(i,j,k,1) = Ajk(1,1);
                        Toep(i,j,k,2) = Ajk(1,2);
                        Toep(i,j,k,3) = Ajk(1,3);
                        Toep(i,j,k,4) = Ajk(2,2);
                        Toep(i,j,k,5) = Ajk(2,3);
                        Toep(i,j,k,6) = Ajk(3,3);
                    end
                    
                end
            else
                rk_to_rj = R1-R0;
                
                rjk = norm(rk_to_rj);
                rjk_hat = (rk_to_rj)/rjk;
                rjkrjk = rjk_hat*rjk_hat';
                
                if abs(rjk)>1e-15
                    Ajk = exp(1i*ko*rjk)/rjk*(ko^2*(I - rjkrjk) + (1i*ko*rjk-1)/rjk^2*(I - 3*rjkrjk)); %Draine & Flatau
                    Toep(i,j,k,1) = Ajk(1,1);
                    Toep(i,j,k,2) = Ajk(1,2);
                    Toep(i,j,k,3) = Ajk(1,3);
                    Toep(i,j,k,4) = Ajk(2,2);
                    Toep(i,j,k,5) = Ajk(2,3);
                    Toep(i,j,k,6) = Ajk(3,3);
                end
            end
            
            
        end
    end
end

opCirculant = circulant_nop_const(Toep);
op_out = fft_operator(opCirculant);



