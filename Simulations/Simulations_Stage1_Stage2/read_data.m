function read_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoCrit = 4;                                 %no. of ITC used in simulations
n = 10;                                     %no. of SP's
%ITC: (1) SBC; (2) FPE; (3) RNML; (4) AICC
vec_N = [600:100:1000 3000 9000 27000];     %sample sizes
Ntr = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%method:
%1 - List (K=5) - This produces Fig. 2
%2 - Inter (K=5) 
%3 - Greedy (K=5)
%4 - List (K=10) - This produces Fig. 3

for method=1:4,
    
    if method<4,
        K = 5;
    else
        K = 10;
    end
    Kbar = K*(K-1)/2;
    
    NS = zeros(length(vec_N),NoCrit);
    ME = zeros(length(vec_N),NoCrit);

    for ind_sp=1:n,
        clear a;
        clear b;
        if method==1,
            fname = strcat('./Results_List_K5/new',num2str(ind_sp));
        elseif method==2,
            fname = strcat('./Results_Inter_K5/new2',num2str(ind_sp));
        elseif method==3,
            fname = strcat('./Results_Greedy_K5/new_3_',num2str(ind_sp));
        else
            fname = strcat('./Results_List_K10/new_4_',num2str(ind_sp));
        end
        load(fname);
        NS = NS + squeeze(sum(a,1));
        ME = ME + squeeze(sum(b,1));
    end
    
    NS = 1 - (NS/Ntr/Kbar/n);
    ME = 1 - (ME/Ntr/Kbar/n);
    
    fig_ind = 2*method-1;
    
    genplot(fig_ind,NS,ME);
end

end%function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function genplot(fig_ind,NS,ME)

AB = zeros(8,8);
AB(:,1:2:end) = NS;
AB(:,2:2:end) = ME;

%Fig. for smaller sample sizes
figure(fig_ind);
bar(1:4,AB(1:4,:),'grouped');
ylim([0 1])
grid on
ax = gca;
set(ax,'XtickLabel',{'600','700','800','900'});
legend('NS-SBC','ME-SBC','NS-FPE','ME-FPE','NS-RNML','ME-RNML','NS-AIC_{c}','ME-AIC_{c}','Location','SouthEast');
xlabel('Sample Size');
ylabel('Similarity index SIM');

%Fig. for larger sample sizes
figure(fig_ind+1);
bar(1:4,AB(5:8,:),'grouped');
ylim([0 1])
grid on
ax = gca;
set(ax,'XtickLabel',{'1000','3000','9000','27000'});
legend('NS-SBC','ME-SBC','NS-FPE','ME-FPE','NS-RNML','ME-RNML','NS-AIC_{c}','ME-AIC_{c}','Location','SouthEast');
xlabel('Sample Size');
ylabel('Similarity index SIM');

end%function
