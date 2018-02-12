function get_plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoCrit = 4;                                 %no. of ITC used in simulations
% ITC = {'SBC', 'FPE', 'RNML', 'AICC'};
AP = {'NO_2','NO','CO','O_3','Rad'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%Init conversion
C = zeros(5,5);
temp = 1;
for i=2:5,
    for j=1:i-1,
        C(i,j) = temp;
        temp = temp+1;
    end
end

%met:
%1 - List
%2 - Mixed List-Greedy
%3 - Greedy

MAT_PHAT = [];
for met=1:3,
    MAT_NS = zeros(11,NoCrit);
    MAT_ME = zeros(11,NoCrit);
    for year=2004:2010,
        results_file = strcat('./Results/Results_Year_',num2str(year),'.mat');   %file to get the results
        load(results_file);
        %Estimated orders
        if met==1,
            MAT_PHAT = [MAT_PHAT; vec_phat];
        end
        %Estimated SP's
        MAT_NS = MAT_NS+get_mat(C,NoCrit,SP_NS,met);
        MAT_ME = MAT_ME+get_mat(C,NoCrit,SP_ME,met);
    end%year
    
    MAT_TOT = zeros(11,2*NoCrit);
    MAT_TOT(:,1:2:end) = MAT_NS;
    MAT_TOT(:,2:2:end) = MAT_ME;
    
    MAT_TOT = MAT_TOT/7; %normalization by the no. of years
    
    MAT_G = [];
    count = 0;
    for i=1:11,
        if sum(MAT_TOT(i,:))>0,
            MAT_G = [MAT_G; MAT_TOT(i,:)];
            if i<11,
                [r,c] = find(C==i);
                label = sprintf('%s%s%s%s%s', '(',AP{r},',',AP{c},')');
                count = count+1;
                labelx{count} = label;
            else
                count = count+1;
                labelx{count} = 'NIL';
            end 
        end
    end
    fig = figure(met);
    colormap(fig,lines);
    bar(1:count,MAT_G,'grouped');
    ylim([0 1])
    grid on
    ax = gca;
    set(ax,'XtickLabel',labelx);
    legend('NS-SBC','ME-SBC','NS-FPE','ME-FPE','NS-RNML','ME-RNML','NS-AIC_{c}','ME-AIC_{c}','Location','NorthWest');
    xlabel('Conditionally independent components');
    ylabel('No. of selections normalized by the no. of years');
end%met

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MAT = get_mat(C,NoCrit,SP,met)

    MAT = zeros(11,NoCrit);
    for crit=1:NoCrit,
        temp = SP{met,crit};
        [r,c,~] = find(temp==0);
        if numel(r)>0,
            ind = 1;
            rind = r(ind);
            cind = c(ind);
            while rind>cind,
                posi = C(rind,cind);
                MAT(posi,crit) = 1;
                ind = ind+1;
                rind = r(ind);
                cind = c(ind);
            end
        else
            MAT(11,crit) = 1;
        end
    end %crit

end% function
