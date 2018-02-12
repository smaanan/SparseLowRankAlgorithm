function mat_time = read_all_results

close all

criteria = {'SBC','FPE','RNML','AICC','Oracle'};
% met = {'ME_3','NS','E_3','S','NS_3','ME','E'};
% met = {'ME_3','NS','E_3','S','NS_3','ME'};
met = {'FL-ME','L-NS','FL-WS','RME-ME','FL-NS','L-ME'};

di{1} = '.\Fast_List\Results_ME3\';
di{2} = '.\Fast_List\Results_NS\';
di{3} = '.\Eichler\Results_Eichler3\';
di{4} = '.\Comparison_Songsiri_JMLR\Results_Sep13\';
di{5} = '.\Fast_List\Results_NS3\';
di{6} = '.\Fast_List\Results_ME\';
% di{7} = '.\Eichler\Results_Eichler\';

N = 1000;
p = 5;

vec_m = [4,10];

mat_time = zeros(5,2);
%Rows: mask4 mask10
%Cols:
%ME3
%NS
%Eichler3
%Songsiri_JMLR
%NS3
%ME
%Eichler

%Results for execution time
for ind_di=1:6,
    for ind_m=1:2,
        ind_mask = vec_m(ind_m);
        mat_time(ind_di,ind_m) = 0;
        for ind_tr=1:9,
            f1 = strcat(di{ind_di},'Res_mask',num2str(ind_mask),'_','p',num2str(p),'_','N',num2str(N),'_','tr',num2str(ind_tr),'.mat');
            FileInfo = dir(f1);
            t1 = datevec(FileInfo.datenum);
            f2 = strcat(di{ind_di},'Res_mask',num2str(ind_mask),'_','p',num2str(p),'_','N',num2str(N),'_','tr',num2str(ind_tr+1),'.mat');
            FileInfo = dir(f2);
            t2 = datevec(FileInfo.datenum);
            mat_time(ind_di,ind_m) = mat_time(ind_di,ind_m) + etime(t2,t1);
        end
        mat_time(ind_di,ind_m) = mat_time(ind_di,ind_m)/9;
    end
end

%Performance
mat_per = read_results(di);

% mat_per{1}

for ii=1:2,
    [time,ind] = sort(mat_time(:,ii),'ascend');
    time = time/10^3;
    perf = mat_per{ii}(ind,:);
%    perf
    figure(ii);
    plot(time,perf(:,1),'rv');
    hold on
    plot(time,perf(:,2),'rp');
    hold on
    plot(time,perf(:,3),'rx');
    hold on
    plot(time,perf(:,4),'ro');
    hold on
    plot(time,perf(:,5),'r^');
    hold on
    
    xlabel('Average runtime for a 10-variate time series (in sec. divided by 10^3)');
    ylabel('Similarity Index SIM');
    legend(criteria);
    
    line([time(1),time(1)],[min(min(perf)),max(max(perf))],'color','g');
    line([time(2),time(2)],[min(min(perf)),max(max(perf))],'color','g');
    line([time(3),time(3)],[min(min(perf)),max(max(perf))],'color','g');
    line([time(4),time(4)],[min(min(perf)),max(max(perf))],'color','g');
    line([time(5),time(5)],[min(min(perf)),max(max(perf))],'color','g');
    line([time(6),time(6)],[min(min(perf)),max(max(perf))],'color','g');
    % line([time(7),time(7)],[min(min(perf)),max(max(perf))],'color','g');
    meti=met(ind);
    
    for t=1:length(time),
        if ii==1,
            if mod(t,2)==0,
                w = 0.6;
            else
                w = 0.55;
            end
        elseif ii==2,
            if mod(t,2)==0,
                w = 0.92;
            else
                w = 0.90;
            end
        else
        end
        text([time(t)-0.03,time(t)-0.03],[w,w],meti{t});
    end
    
   
end %ii


end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mat_per = read_results(di)

NoCrit = 4;

vec_N = 1000;
Nmax   = max(vec_N);

%model orders; in estimation, we use the true order
vec_p = 5;
%no. of variables:
K      = 10;
kbar = K*(K-1)/2;
%no. of trials:
Ntr    = 10;



for ind_di=1:size(di,2),

ii = 0;
for ind_mask=[4,10],
    ii = ii+1;
    Mat{ii} = zeros(1,5);
    for ind_p = 1:length(vec_p), %one single value for p
        p = vec_p(ind_p);
        for ind_tr = 1:Ntr,
            for ind_N = 1:length(vec_N), %one single value for N
                N = vec_N(ind_N);
                fres = strcat(di{ind_di},'Res_mask',num2str(ind_mask),'_','p',num2str(p),'_','N',num2str(N),'_','tr',num2str(ind_tr),'.mat');
                load(fres);
                
                if ind_di == 1,
                    vec_crit = vec_crit_me;
                    vec_mask = vec_mask_me;
                elseif ind_di == 2,
                    vec_crit = vec_crit_ns;
                    vec_mask = vec_mask_ns;
                elseif ind_di == 3,
                    vec_crit = vec_crit_ei;
                    vec_mask = vec_mask_ei;
                elseif ind_di == 4,
                elseif ind_di == 5,
                    vec_crit = vec_crit_ns;
                    vec_mask = vec_mask_ns;
                elseif ind_di == 6,
                    vec_crit = vec_crit_me;
                    vec_mask = vec_mask_me;
                elseif ind_di == 7,
                    vec_crit = vec_crit_ei;
                    vec_mask = vec_mask_ei;
                end
                        
                for cc=1:4,
                    [~,w] = min(vec_crit(cc,:));
                    Mat{ii}(ind_p,cc) = Mat{ii}(ind_p,cc)+sum(sum(xor(mask_true,vec_mask{w})))/2;
                end %cc
                oracle = Inf;
                for jj=1:size(vec_mask,2);
                    mask = vec_mask{jj};
                    temp = sum(sum(xor(mask,mask_true)))/2;
                    oracle = min(oracle,temp);
                end
                Mat{ii}(ind_p,5) = oracle;
            end %N
        end %tr
    end %p
mat_per{ii}(ind_di,:) = ones(size(Mat{ii})) - (Mat{ii}/kbar/Ntr);
end %mask

end %ind_di


        
end %function



