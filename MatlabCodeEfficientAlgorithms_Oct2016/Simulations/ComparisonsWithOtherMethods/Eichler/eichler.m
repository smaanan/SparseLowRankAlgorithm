function [A, Q] = eichler(cur_mask,mask_estim,Rest,p,data,Noit,Ngrid,rat0)

[~,K] = size(data);
vec_om = 0:(pi/Ngrid):pi;
noz = sum(sum(mask_estim))/2;
kbar = K*(K-1)/2;
Resti = Rest;
it = 0;

% Nmax = 100000;
% N = 27000;
% N = 1000;
% K = 10;
% p = 5;
% Ngrid = 1024;
% ind_mask = 4;
% load('mask_set.mat'); 
% Noit = 50;
%%%%%%%%%

% mask = mask_set(:,(ind_mask-1)*K+1:ind_mask*K);
% cur_mask = mask;
% [~, A, Sigma] = gen_IS_model(K, p, mask);
% all_data = arsim(zeros(K,1), A, Sigma, Nmax);
% save('temp.mat','A','Sigma','all_data','cur_mask');

% load('temp.mat');
% load('Data_mask10_p5_tr10.mat');


       
% data = all_data(1:N,:);
% Rest = comp_R(data, p);


% vec_sbc = zeros(1,Noit);
% vec_sumz = zeros(1,Noit);
% vec_sumnz = zeros(1,Noit);




while it<Noit,  
    it = it+1;
   
    [A,~] = solve_yw(Resti, K, p);
    % vec_sbc(it) = comp_sbc(data, Sigma, p, noz);
    
    Q = comp_Q(A, K, p);
    
    sumz = 0;
    sumnz = 0;
    for u=1:p+1,
        sumz = sumz + sum(sum(triu(abs(Q{u}.*mask_estim),1)))/noz/(p+1);
%         if (kbar-noz)>0,
%             sumnz = sumnz + sum(sum(triu(abs(Q{u}.*cur_mask),1)))/(kbar-noz)/(p+1);
%         end
    end
%     rat = sumnz/sumz;
    if ( (it>3) && (sumz<10^(-4)) ),
%     if ( (it>2) && ((rat<rat0) && (sumnz>0)) || ((sumz<10^(-3)) && (sumnz==0)) ),
        break;
    end
%     vec_sumz(it) = sumz;
%     vec_sumnz(it) = sumnz;
  
    [S, Sinv] = comp_S_Sinv(Q, vec_om);
    
    for ii=1:K,
        for jj=ii+1:K,
            if cur_mask(ii,jj)==0,
                for ind_om = 1:length(vec_om),
                    [S{ind_om}, Sinv{ind_om}] = invest(S{ind_om}, Sinv{ind_om}, ii, jj, K);
                end
            end
        end %jj
    end %ii

    Resti = comp_Rest_int(S, vec_om, Ngrid, p);
    
end %it

end %function



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, Sinv_tot] = comp_S_Sinv(Yestim, vec_om)

i = sqrt(-1);

Sinv0 = Yestim{1};
Sinv = Sinv0;
for ind=1:length(vec_om),
    om = vec_om(ind);
    for ell=2:size(Yestim,2),
        Sinv = Sinv + exp(-i*om*(ell-1))*Yestim{ell} +exp(i*om*(ell-1))*transpose(Yestim{ell});
    end
    Sinv_tot{ind} = Sinv;
    S{ind} = inv(Sinv);
    Sinv = Sinv0;
end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [M, Minv] = invest(S, Sinv, ii, jj, K)

% jj>ii

D = Sinv(ii,ii)*Sinv(jj,jj) - Sinv(ii,jj)^2;

M = S;
M(ii,jj) = S(ii,jj) + Sinv(ii,jj)/D;
M(jj,ii) = conj(M(ii,jj));

Minv = zeros(K,K);
for r=1:K,
    for c=r:K,
        if r==ii,
            if c==jj,
                Minv(r,c) = 0;
            elseif c==ii,
                Minv(r,c) = D/Sinv(jj,jj);
            else
                Minv(r,c) = Sinv(ii,c) - Sinv(ii,jj)*Sinv(jj,c)/Sinv(jj,jj);
            end
        elseif r==jj,
            if c==jj,
                Minv(r,c) = D/Sinv(ii,ii);
            else
                Minv(r,c) = Sinv(jj,c) - Sinv(ii,jj)*Sinv(ii,c)/Sinv(ii,ii);
            end
        else
            t1 = Sinv(jj,c)-Sinv(ii,jj)*Sinv(ii,c)/Sinv(ii,ii);
            t2 = Sinv(ii,c)-Sinv(ii,jj)*Sinv(jj,c)/Sinv(jj,jj);
            Minv(r,c) = Sinv(r,c) - (Sinv(ii,jj)/D)*( Sinv(ii,r)*t1 + Sinv(jj,r)*t2 );
        end
    end %c
end %r 

for r=1:K,
    Minv(r,r) = real(Minv(r,r));
end

Minv = Minv + conj(transpose(triu(Minv,1)));

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Resti = comp_Rest_int(S, vec_om, Ngrid, p)

i = sqrt(-1);

for u=0:p,
    temp = 0;
    for ind_om = 1:(length(vec_om)-1),
        om = vec_om(ind_om);
        temp = temp+2*real(S{ind_om}*exp(i*om*u));
    end
    Resti{u+1} = temp*(pi/Ngrid)/(2*pi);
end

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function sbc = comp_sbc(data, Sigma, phat, noz)
% 
% N = size(data, 1);
% T = N - phat;
% n = size(data,2);
% 
% %No. of param. [Songsiri et al., book chapter, p. 106]
% nop = n*(n+1)/2 - noz + phat*(n^2-2*noz);
% nop1= nop/n;
%  
% %SBC
% sbc = T*log(det(Sigma)) + nop*log(T);
% 
% end %function
