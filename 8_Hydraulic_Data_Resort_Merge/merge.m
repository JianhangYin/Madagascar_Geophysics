clear all
close all
addpath('D:\Paper\Data\CottonValley\')
cd('resorted');
daysinfo = dir('*');
for days = 5:1:(length(daysinfo))
    cd(daysinfo(days).name);  
    compinfo = dir('*');
    mkdir(daysinfo(days).name);
    for comp = 3:1:(length(compinfo))
        cd(compinfo(comp).name);
        statinfo = dir('*');
        for stat = 3:1:(length(statinfo))
            disp(['Working on day: ' daysinfo(days).name ', component: ' compinfo(comp).name ', station: ' statinfo(stat).name])
            cd(statinfo(stat).name);
            tic
            dirinfo = dir('*.sac');
            data_files = cell(1,1);
            count = 0;
            merged = zeros(86400000,1);
            coeffi = zeros(86400000,1);
            for ifile = 1:1:(length(dirinfo))
                if (dirinfo(ifile).bytes < 32)
                    disp([dirinfo(ifile).name 'is demaged']);
                    continue;
                end
                fname=dirinfo(ifile).name;
                data=rsac('big-endian',fname);
                count=count+1;
                data_files{count}=data;
                kztime = sum(data_files{count}(73:76,3) .* [3600,60,1,0.001]');
                data_files{count}(:,1) = data_files{count}(:,1) + kztime;
                begintime = round(data_files{count}(1,1)*1000);
                numbetime = length(data_files{count}(:,1));
                merged(begintime:begintime+numbetime-1) = data_files{count}(:,2).*(1-coeffi(begintime:begintime+numbetime-1))+merged(begintime:begintime+numbetime-1).*coeffi(begintime:begintime+numbetime-1);
                coeffi(begintime:begintime+numbetime-1) = 0.5;
            end
            cd('..');
            %save(['..\' daysinfo(days).name '\' compinfo(comp).name '_' statinfo(stat).name '.mat'], 'merged');
            fid = fopen(['..\' daysinfo(days).name '\' compinfo(comp).name '_' statinfo(stat).name '.rsf@'],'w');
            fwrite(fid,merged,'float');
            toc
        end
        cd('..');
    end
    cd('..')
end











