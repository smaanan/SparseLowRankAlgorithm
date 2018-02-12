function get_file_time


% Get individual components of date & time in 1 Sec resolution
FileInfo = dir('Res_mask4_p5_N1000_tr2.mat');
[Y, M, D, H, MN, S] = datevec(FileInfo.datenum);

FileInfo = dir('Res_mask4_p5_N1000_tr1.mat');
t1 = datevec(FileInfo.datenum);

FileInfo = dir('Res_mask4_p5_N1000_tr2.mat');
t2 = datevec(FileInfo.datenum);

etime(t2,t1)

end

