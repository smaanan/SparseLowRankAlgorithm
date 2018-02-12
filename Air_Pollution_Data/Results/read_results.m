function read_results( year )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoCrit = 4;                                 %no. of ITC used in simulations
ITC = {'SBC', 'FPE', 'RNML', 'AICC'};
NoMet = 3;                                  %no. of methods
MET = {'List','L-G','Greedy'};
AP = {'NO$_2$','NO','CO','O$_3$','Rad'};

r1 = '\begin{table}[h!]';
r2 = '\begin{center}';
r3 = '\begin{tabularx}{1.0\textwidth}{>{\hsize=0.2\hsize}X>{\centering\hsize=0.2\hsize}X>{\centering\hsize=0.2\hsize}X>{\centering\hsize=0.2\hsize}X>{\hsize=0.2\hsize}X}';
r4 = '\hline';
r5 = '\multicolumn{5}{l}{Estimated orders for the VAR model  }\\'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

results_file = strcat('Results_Year_',num2str(year),'.mat');   %file to get the results
load(results_file);

results_file_out = strcat('Results_Year_',num2str(year),'.tex');   %file to save the results

fid = fopen(results_file_out,'w');

fprintf(fid,'%s \n',r1);
fprintf(fid,'%s \n',r2);
fprintf(fid,'%s \n',r3);
fprintf(fid,'%s \n',r4);
fprintf(fid,'%s \n',r5);

%Estimated orders
fprintf(fid,'%s','&');
for crit=1:NoCrit-1,
    fprintf(fid,'%s %s',ITC{crit},'&');
end
fprintf(fid,'%s %s \n',ITC{NoCrit},'\\');
fprintf(fid,'%s','&');
for crit=1:NoCrit-1,
     fprintf(fid,'%2i %s',vec_phat(crit),'&');
end
fprintf(fid,'%2i %s \n',vec_phat(NoCrit),'\\');    

%Estimated SP's
for opt=1:2,
    if opt==1,
        fprintf(fid,'\n %s \n %s \n', '\hline ', '\multicolumn{5}{l}{Near Sparse}\\');
        SP = SP_NS;
        write_text(fid,NoCrit,ITC,NoMet,MET,AP,SP);
    elseif opt==2,
        fprintf(fid,'\n %s \n %s \n', '\hline ', '\multicolumn{5}{l}{Maximum Entropy}\\');
        SP = SP_ME;
        write_text(fid,NoCrit,ITC,NoMet,MET,AP,SP);
    else
    end
end%case

fprintf(fid,'\n %s \n %s \n %s \n', '\hline', '\end{tabularx}', '\end{center}');
fprintf(fid,'%s %s %s \n %s \n ', '\caption{Year', num2str(year),  '}', '\end{table}'); 

fclose(fid);

end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function write_text(fid,NoCrit,ITC,NoMet,MET,AP,SP)

fprintf(fid,'%s','&');
for crit=1:NoCrit-1,
    fprintf(fid,'%s %s',ITC{crit},'&');
end
fprintf(fid,'%s %s \n',ITC{NoCrit},'\\');

for met=1:NoMet,
    fprintf(fid,'%s %s',MET{met},'&');
    for crit=1:NoCrit,
        temp = SP{met,crit};
        [r,c,~] = find(temp==0);
        if numel(r)>0,
            ind = 1;
            rind = r(ind);
            cind = c(ind);
            while rind>cind,
                fprintf(fid,'%s%s%s%s%s',' (',AP{rind},',',AP{cind},')');
                ind = ind+1;
                rind = r(ind);
                cind = c(ind);
            end
        else
            fprintf(fid,'%s','None');
        end
        if crit<NoCrit,
            fprintf(fid,'%s','&');
        elseif crit==NoCrit,
            fprintf(fid,'%s \n','\\');
        end 
    end %crit
end %met


end% function
