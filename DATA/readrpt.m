if isempty(legend)
    clear
end
clc
R=1;
N=1;
w=[0.01,5000];
[FileName,PathName] = uigetfile('*.rpt','请选择报告文件：');
Name=FileName(1:end-4);
h=0;
if ~isempty(legend)
    for i=1:length(lstr)
        if length(Name)==length(char(lstr(i))) && strcmpi(Name,char(lstr(i)))
            h=1;
            break
        end
    end
end
clc
if h==0
    fid=fopen([PathName,FileName],'rt+');
    for i=1:32
        h=fgets(fid);
        h(h(1:min(32,length(h)))==' ')=[];
        if h(1)=='X'
            break
        end
    end
    data=textscan(fid,' %f %f \n','HeaderLines',2);
    fclose(fid);
    data=[cell2mat(data(1)) cell2mat(data(2)) 20*log10(cell2mat(data(2)))/N];
    h=data(:,3);
    h(data(:,2)<1)=NaN;
    data=[data h];
    %eval([Name,'=data;']);
    hold on
    plot(data(data(:,1)>=w(1)&data(:,1)<=w(2),1),data(data(:,1)>=w(1)&data(:,1)<=w(2),4),'LineWidth',3)
    if isempty(legend)
        lstr={Name};
        xlabel('频率/Hz')
        ylabel('响应/dB')
        title('频响曲线：')
        grid on
    else lstr=[lstr cellstr(Name)];
    end
    R=cell(1,length(lstr));
    for i=1:length(lstr)
        h=char(lstr(i));
        h(h=='_')='-';
        R(i)=cellstr(h);
    end
    legend(char(R))
    hold off
    clear ans data fid
end
clear i Name FileName PathName R w h