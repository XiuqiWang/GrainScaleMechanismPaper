close all
clear

% 指定文件名的基本部分和扩展名
filerootname = 'S005M10LBIni';
fileExt = '.vtu';
%%you need to change these each time!!!
N=2725;
Nstart=204;
Nend=704;

fileBase = ['E:\cluster\WilletViscous\LBIniTP\S005\',filerootname,'Particle_'];

Position = {};
Velocity = {};
Radii = {};
FullLiquidVolume = {};
LiquidFilmVolume = {};
LiquidBridgeVolume = {};

% 循环读取文件
for i = Nstart:Nend  % 假设有10个文件，文件名为 data_1.txt 到 data_10.txt
    dataCell = {};
    % 构建当前文件名
    filename = [fileBase, num2str(i), fileExt];
    
    filein = filename;
    fid = fopen(filein,'r');
    
    while ~feof(fid)%till the end of a file
    tline=fgetl(fid);
    if contains(tline,'<')%contains '<'
        continue
    else
        dataCell{end+1} = tline;
        %fprintf(fout,'%s\n',tline);
    end
    end
    fclose(fid);
    %fclose(fout);
    
    Position{i-Nstart+1}=dataCell(1:N);
    Velocity{i-Nstart+1}=dataCell(N+1:N*2);
    Radii{i-Nstart+1}=dataCell(N*3+1:N*4);
    FullLiquidVolume{i-Nstart+1}=dataCell(N*5+1:N*6);
    LiquidFilmVolume{i-Nstart+1}=dataCell(N*6+1:N*7);
    LiquidBridgeVolume{i-Nstart+1}=dataCell(N*7+1:N*8);
end

save([filerootname,'.mat'],'Position','Velocity','Radii','FullLiquidVolume','LiquidFilmVolume','LiquidBridgeVolume');
%save([filerootname,'.mat'],'Position','Velocity','Radii');
