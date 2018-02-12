function [totalp1, totalp2] = draw_graphs_order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots the results produced by run_exp_order.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
load('Result_Order_Estimation');

half = length(vec_N)/2;
vec_N1 = vec_N(1:half);
ord_eq1 = ord_eq(:,1:half);
[totalp1, ~] = comp_total(ord_eq1, vec_p, vec_N1, Ntr);
vec_N2 = vec_N(half+1:end);
ord_eq2 = ord_eq(:,half+1:end);
[totalp2, ~] = comp_total(ord_eq2, vec_p, vec_N2, Ntr);

fig_no = 1;
for ind_p=1:length(vec_p),
    figure(fig_no);
    fig_no = fig_no+1;
    
    bar(vec_N1,...
        totalp1(:,(ind_p-1)*length(vec_N1)+1:ind_p*length(vec_N1))',...
        'grouped');
    xlabel('Sample size');
    ylabel('Empirical prob. of correctly selecting the order');
    legend('SBC','FPE','RNML','AIC','AICc','KIC','KICc','Location', 'Best');
    tit = strcat('p^\circ=', num2str(vec_p(ind_p)));
    title(tit);
    
    figure(fig_no);
    fig_no = fig_no+1;
    bar(log(vec_N2), ...
        totalp2(:,(ind_p-1)*length(vec_N2)+1:ind_p*length(vec_N2))', ...
        'grouped');
    xlabel('Sample size (logarithmic scale)');
    labels = cell(1,4);
    for ii=1:half,
        labels{ii} = num2str(vec_N(half+ii));
    end
    set(gca,'XTick',log(vec_N2));
    set(gca,'XTickLabel',labels);
    ylabel('Empirical prob. of correctly selecting the order');
    legend('SBC','FPE','RNML','AIC','AICc','KIC','KICc','Location', 'Best');
    tit = strcat('p^\circ=', num2str(vec_p(ind_p)));
    title(tit);
end
end  %function draw_graphs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [totalp, total] = comp_total(ord_eq, vec_p, vec_N, Ntr)

%Graphs total = for order selection
[~, Nmask] = size(ord_eq{1,1});

sbc = zeros(length(vec_p),length(vec_N));
fpe = zeros(length(vec_p),length(vec_N));
rnml = zeros(length(vec_p),length(vec_N));
aic = zeros(length(vec_p),length(vec_N));
aicc = zeros(length(vec_p),length(vec_N));
kic = zeros(length(vec_p),length(vec_N));
kicc = zeros(length(vec_p),length(vec_N));

for ind_N=1:length(vec_N),
    temp1 = ord_eq{1,ind_N};
    temp2 = ord_eq{2,ind_N};
    temp3 = ord_eq{3,ind_N};
    temp4 = ord_eq{4,ind_N};
    temp5 = ord_eq{5,ind_N};
    temp6 = ord_eq{6,ind_N};
    temp7 = ord_eq{7,ind_N};
    
    %Sum over all masks
    for ind_p=1:length(vec_p),
        sbc(ind_p,ind_N) = sum(temp1(ind_p,:))/Nmask/Ntr;
        fpe(ind_p,ind_N) = sum(temp2(ind_p,:))/Nmask/Ntr;
        rnml(ind_p,ind_N) = sum(temp3(ind_p,:))/Nmask/Ntr;
        aic(ind_p,ind_N) = sum(temp4(ind_p,:))/Nmask/Ntr;
        aicc(ind_p,ind_N) = sum(temp5(ind_p,:))/Nmask/Ntr;
        kic(ind_p,ind_N) = sum(temp6(ind_p,:))/Nmask/Ntr;
        kicc(ind_p,ind_N) = sum(temp7(ind_p,:))/Nmask/Ntr;
    end
end

sbct = sbc';
fpet = fpe';
rnmlt = rnml';
aict = aic';
aicct = aicc';
kict = kic';
kicct = kicc';

totalp = [sbct(:) fpet(:) rnmlt(:) aict(:) aicct(:) kict(:) kicct(:)]';


%Sum over all "true" orders
total = [
    sum(sbc,1)/size(sbc,1)
    sum(fpe,1)/size(fpe,1)
    sum(rnml,1)/size(rnml,1)
    sum(aic,1)/size(aic,1)
    sum(aicc,1)/size(aicc,1)
    sum(kic,1)/size(kic,1)
    sum(kicc,1)/size(kicc,1)];
end %function comp_total