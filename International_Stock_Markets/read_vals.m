function read_vals

close all

ITC = {'SBC','RNML'};

plot_vals(ITC{1});
plot_vals(ITC{2});

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_vals(crit)
%Results from Abdelwahab et al.
[true_mat_weak, true_mat_strong] = build_true_mat;
k = size(true_mat_weak,1);

%Results from experiments
fname = strcat('Results_FinData_Vals_Masks_',crit,'.mat');
load(fname);

valse = [val0; vals];
[mm, ww] = min(valse);
valse = valse-mm;

sred = [];
cred = [];
tgr = [];

for i=2:length(valse),
    if i==2,
        temp = xor(masks{i-1,1},ones(k,k));
    else
        temp = xor(masks{i-1,1},masks{i-2,1});
    end
    
    if sum(sum(and(temp,true_mat_weak)))==2,
        if sum(sum(and(temp,true_mat_strong)))==2,
            sred = [sred; i,valse(i)];
        else
            cred = [cred; i,valse(i)];
        end
    else
        tgr = [tgr; i,valse(i)];
    end            
end %for

figure;
plot([1:length(valse)]-1, valse, '-k', ...
    sred(:,1)-1, sred(:,2), '*r', ...
    cred(:,1)-1, cred(:,2), 'ro', ...
    tgr(:,1)-1, tgr(:,2), 'vg');
legend('ITC value','strong PSC','weak PSC','zero PSC','Location','NorthWest');
ylabel(crit);
xlabel('Index (i) in List procedure');
hold on
line([ww,ww]-1,[min(valse),max(valse)],'LineStyle','--','Color','k');
ylim([0 max(valse)]);

end %function

