%clearvars -except FifteenMinHeLa FifteenMinRPE1 FifteenMinWTDose FifteenMinN386CDose FifteenMinshSPAST
clear all
close all

% setTx=0; %set to 0 if DMSO, 1 is spastazoline treated

testTitle = '20201130_HeLa';

% Specify the folder where the files live.
DtheFolder = 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_SPASTMT\Tracking\Image_files\Tracks\MATLAB_Data\HeLa\DMSO';
StheFolder = 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_SPASTMT\Tracking\Image_files\Tracks\MATLAB_Data\HeLa\SPZ';

%Location files to be saved
cd('D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_SPASTMT\Tracking\Image_files\Tracks\MATLAB_Data\Analysis');



% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(DtheFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', DtheFolder);
  uiwait(warndlg(errorMessage));
  return;
end

if ~isdir(StheFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', StheFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name
% pattern.c
DfPattern = fullfile(DtheFolder, '*.mat'); % Change to whatever pattern you need.
DmatFiles = dir(DfPattern);

SfPattern = fullfile(StheFolder, '*.mat'); % Change to whatever pattern you need.
SmatFiles = dir(SfPattern);

AllData=[];
AllDMSO=[];
AllSPZ=[];

AllKin=[];
AllMid=[];

AllDKin=[];
AllSKin=[];
AllDMid=[];
AllSMid=[];

AllTub=[];
AllNoTub=[];

AllDTub=[];
AllDNoTub=[];
AllSTub=[];
AllSNoTub=[];

d10spz=['10 ',char(181),'M spastazoline'];
micron=[char(181),'m'];

for n = 1:length(DmatFiles) 
    clear KinCore MidzoneCore Noncore KinCore_Tub MidzoneCore_Tub Noncore_Tub KinCore_noTub MidzoneCore_noTub Noncore_noTub
    
    baseName = DmatFiles(n).name;
    fullName = fullfile(DtheFolder, baseName);
    [filepath,name,ext] = fileparts(fullName);
    matName=name;
    
    load(fullfile(filepath, strcat(name,ext)));
    
    if exist('KinCore','var')==1
        AllDMSO=[AllDMSO;KinCore];
        AllDKin=[AllDKin;KinCore];
        if exist('KinCore_Tub','var')==1
            AllDTub=[AllDTub;KinCore_Tub];
        else
        end
        if exist('KinCore_noTub','var')==1
            AllDNoTub=[AllDNoTub;KinCore_noTub];
        else
        end
    else
    end
    
    if exist('MidzoneCore','var')==1
        AllDMSO=[AllDMSO;MidzoneCore];
        AllDMid=[AllDMid;MidzoneCore];
        if exist('MidzoneCore_Tub','var')==1
            AllDTub=[AllDTub;MidzoneCore_Tub];
        else
        end
        if exist('MidzoneCore_noTub','var')==1
            AllDNoTub=[AllDNoTub;MidzoneCore_noTub];
        else
        end
    else
    end
    
    if exist('Noncore','var')==1
        AllDMSO=[AllDMSO;Noncore];
        if exist('Noncore_Tub','var')==1
            AllDTub=[AllDTub;Noncore_Tub];
        else
        end
        if exist('Noncore_noTub','var')==1
            AllDNoTub=[AllDNoTub;Noncore_noTub];
        else
        end
    else
    end
    
end

for n = 1:length(SmatFiles) 
    clear KinCore MidzoneCore Noncore KinCore_Tub MidzoneCore_Tub Noncore_Tub KinCore_noTub MidzoneCore_noTub Noncore_noTub
    
    baseName = SmatFiles(n).name;
    fullName = fullfile(StheFolder, baseName);
    [filepath,name,ext] = fileparts(fullName);
    matName=name;
    
    load(fullfile(filepath, strcat(name,ext)));

    if exist('KinCore','var')==1
        AllSPZ=[AllSPZ;KinCore];
        AllSKin=[AllSKin;KinCore];
        if exist('KinCore_Tub','var')==1
            AllSTub=[AllSTub;KinCore_Tub];
        else
        end
        if exist('KinCore_noTub','var')==1
            AllSNoTub=[AllSNoTub;KinCore_noTub];
        else
        end
    else
    end
    
    if exist('MidzoneCore','var')==1
        AllSPZ=[AllSPZ;MidzoneCore];
        AllSMid=[AllSMid;MidzoneCore];
        if exist('MidzoneCore_Tub','var')==1
            AllSTub=[AllSTub;MidzoneCore_Tub];
        else
        end
        if exist('MidzoneCore_noTub','var')==1
            AllSNoTub=[AllSNoTub;MidzoneCore_noTub];
        else
        end
    else
    end
    
    if exist('Noncore','var')==1
        AllSPZ=[AllSPZ;Noncore];
        if exist('Noncore_Tub','var')==1
            AllSTub=[AllSTub;Noncore_Tub];
        else
        end
        if exist('Noncore_noTub','var')==1
            AllSNoTub=[AllSNoTub;Noncore_noTub];
        else
        end
    else
    end
    
end
    
AllData=[AllData;AllDMSO];
AllData=[AllData;AllSPZ];

AllKin=[AllKin;AllDKin];
AllKin=[AllKin;AllSKin];

AllMid=[AllMid;AllDMid];
AllMid=[AllMid;AllSMid];

AllTub=[AllTub;AllDTub];
AllTub=[AllTub;AllSTub];

AllNoTub=[AllNoTub;AllDNoTub];
AllNoTub=[AllNoTub;AllSNoTub];

save(testTitle);

%(1) mean spas signal, (2) max tub signal, (3) spas track lifetime, (4) trackmate ID
%%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%

h1=figure;
title('Tubulin Intensity vs. Lifetime')
hold on
    x=double(AllData(:,3));
    y=double(AllData(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);

    
h2=figure;
title('Spastin Intensity vs. Lifetime')
hold on
    x=double(AllData(:,3));
    y=double(AllData(:,1));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('mean spastin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
    
h3=figure;
title('Tubulin Intensity vs. Spastin Intensity')
hold on
    x=double(AllData(:,1));
    y=double(AllData(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('mean spastin intensity (a.u.)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
%%%%%%%%%%%%%%%%%%%DMSO%%%%%%%%%%%%%%%%%%

h4=figure;
title('Tubulin Intensity vs. Lifetime')
hold on
    x=double(AllDMSO(:,3));
    y=double(AllDMSO(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);

    
h5=figure;
title('Spastin Intensity vs. Lifetime')
hold on
    x=double(AllDMSO(:,3));
    y=double(AllDMSO(:,1));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('mean spastin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
    
h6=figure;
title('Tubulin Intensity vs. Spastin Intensity')
hold on
    x=double(AllDMSO(:,1));
    y=double(AllDMSO(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('mean spastin intensity (a.u.)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
    %%%%%%%%%%%%%SPZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
h7=figure;
title('Tubulin Intensity vs. Lifetime')
hold on
    x=double(AllSPZ(:,3));
    y=double(AllSPZ(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','b','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);

    
h8=figure;
title('Spastin Intensity vs. Lifetime')
hold on
    x=double(AllSPZ(:,3));
    y=double(AllSPZ(:,1));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','b','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('mean spastin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
    
h9=figure;
title('Tubulin Intensity vs. Spastin Intensity')
hold on
    x=double(AllSPZ(:,1));
    y=double(AllSPZ(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','b','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('mean spastin intensity (a.u.)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
    %%%%%%%%%%%%%%%%%%%TUB ONLY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
h10=figure;
title('Tubulin Intensity vs. Lifetime')
hold on
    x=double(AllDTub(:,3));
    y=double(AllDTub(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);

    
h11=figure;
title('Spastin Intensity vs. Lifetime')
hold on
    x=double(AllDTub(:,3));
    y=double(AllDTub(:,1));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('spastin track lifetime (frames)') 
    ylabel(strcat('mean spastin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
    
h12=figure;
title('Tubulin Intensity vs. Spastin Intensity')
hold on
    x=double(AllDTub(:,1));
    y=double(AllDTub(:,2));
    mdl=fitlm(x,y);
    Rord=num2str(round(mdl.Rsquared.Ordinary,2));
    scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
    %plot(mdl)
    xlabel('mean spastin intensity (a.u.)') 
    ylabel(strcat('maximum tubulin intensity (a.u.)'))
    Rtext=strcat('R^{2}=',Rord);
    annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
    set(gca,'FontSize', 16);
    
    %%%%%%%%%%%%%%%%%NO TUB%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
%     h13=figure;
% title('Tubulin Intensity vs. Lifetime')
% hold on
%     x=double(AllDNoTub(:,3));
%     y=double(AllDNoTub(:,2));
%     mdl=fitlm(x,y);
%     Rord=num2str(round(mdl.Rsquared.Ordinary,2));
%     scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
%     %plot(mdl)
%     xlabel('spastin track lifetime (frames)') 
%     ylabel(strcat('maximum tubulin intensity (a.u.)'))
%     Rtext=strcat('R^{2}=',Rord);
%     annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
%     set(gca,'FontSize', 16);
% 
%     
% h14=figure;
% title('Spastin Intensity vs. Lifetime')
% hold on
%     x=double(AllDNoTub(:,3));
%     y=double(AllDNoTub(:,1));
%     mdl=fitlm(x,y);
%     Rord=num2str(round(mdl.Rsquared.Ordinary,2));
%     scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
%     %plot(mdl)
%     xlabel('spastin track lifetime (frames)') 
%     ylabel(strcat('mean spastin intensity (a.u.)'))
%     Rtext=strcat('R^{2}=',Rord);
%     annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
%     set(gca,'FontSize', 16);
%     
%     
% h15=figure;
% title('Tubulin Intensity vs. Spastin Intensity')
% hold on
%     x=double(AllDNoTub(:,1));
%     y=double(AllDNoTub(:,2));
%     mdl=fitlm(x,y);
%     Rord=num2str(round(mdl.Rsquared.Ordinary,2));
%     scatter(x,y,'filled','k','MarkerFaceAlpha',0.2)
%     %plot(mdl)
%     xlabel('mean spastin intensity (a.u.)') 
%     ylabel(strcat('maximum tubulin intensity (a.u.)'))
%     Rtext=strcat('R^{2}=',Rord);
%     annotation('textbox',[0.72 0.42 0.5 0.5],'Linestyle','none','String',Rtext,'Fontsize',16)
%     set(gca,'FontSize', 16);
%     