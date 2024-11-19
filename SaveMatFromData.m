close all
% clear
% data0 = read_data('S004DryLBIni.data',14,[1 6]);
% data1 = read_data('S004M1LBR1Ini.data',14,[2 7]);
% data2 = read_data('S004M5LBIni.data',14,[2 7]);
% data3 = read_data('S004M10LBIni.data',14,[2 7]);
% data4 = read_data('S004M20LBIni.data',14,[2 7]);
% data5 = read_data('FTSS005M20LBIni.data',14,[2 7]);
% data6 = read_data('FiftyTSS005M20LBIni.data',14,[2 7]);

%save data to .mat for future use
save('FiftyTSS005M20LBIni.mat','data6', '-v7.3');
