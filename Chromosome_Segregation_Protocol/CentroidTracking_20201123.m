clear all; close all;
tic
mincheck=[];
iii=1;
% Specify the folder where the files live.
myFolder = 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_Chromosome_Separation\CellLines_wekaAnalysis\RPE1Parental\SPZ\Analysis\DNAresults';
analysisTitle = '2_20210216_RPE1_SPZ';
minFr=9;
maxFr=50;

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name
% pattern.c
filePattern = fullfile(myFolder, '*.csv'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

for k = 1 : length(theFiles)
  close all;
  clearvars -except myFolder filePattern theFiles k allDist analysisTitle mincheck iii minFr maxFr; 
  baseFileName = theFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  [filepath,name,ext] = fileparts(fullFileName);
  splitName = split(name,'_uncorr');
  %ROIext = '.roi';
  %strFilename = (fullfile(filepath, strcat(splitName{1,1},'_ROI',ROIext)));

    Data=readtable(fullfile(filepath, strcat(name,ext)));
    %ROI = ReadImageJROI(strFilename);
    ROIs=[table2array(Data(:,14)), table2array(Data(:,4)), table2array(Data(:,5)), table2array(Data(:,3)), table2array(Data(:,8)), table2array(Data(:,13))]; %frame,x,y,area,perimeter,circularity
    cd('D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_Chromosome_Separation\CellLines_wekaAnalysis\MatLabAnalyses')

    %filter distances
    int=1; %interval in (min)
    pxl = 1;
    %pxl = (1/9.0909); %converts pixels to microns
    one=[];two=[];three=[];four=[];five=[];six=[];

    for i=0:length(ROIs(:,1))+1

        target=find(ROIs(:,1)==i);

       if length(target)==2 
        one=cat(1,one,ROIs(target,1));
        two=cat(1,two,ROIs(target,2));
        three=cat(1,three,ROIs(target,3));
        four=cat(1,four,ROIs(target,4));
        five=cat(1,five,ROIs(target,5));
        six=cat(1,six,ROIs(target,6));
        %% one, two three and four are filtered
       else
       clear target 
       end
    end

    % centroidDistFilt are the filtered distances
    centroidDistFilt(:,1)=one*int; %time (min)
    centroidDistFilt(:,2)=two; %x (um)
    centroidDistFilt(:,3)=three; %y (um)
    centroidDistFilt(:,4)=four; %area (um)
    centroidDistFilt(:,5)=five; %perimeter (um)
    centroidDistFilt(:,6)=six; %circularity

    clear one two three
    count=1;

    for i=1:length(centroidDistFilt(:,1))
        time=centroidDistFilt(i,1);
        target=find(centroidDistFilt(:,1)==time);

        if length(target)==2
        distances{count,1}(:,1)=centroidDistFilt(target,1); % time (min)
        distances{count,1}(:,2)=centroidDistFilt(target,2); % x in mum
        distances{count,1}(:,3)=centroidDistFilt(target,3); % y in mum
        distances{count,1}(:,4)=centroidDistFilt(target,4); % area (um)
        distances{count,1}(:,5)=centroidDistFilt(target,5); % perimeter (um)
        distances{count,1}(:,6)=centroidDistFilt(target,6); % circularity
        count=count+1;
        else
        end

    end

    calcDist = [];
    j=1;
    for i=1:length(distances)
        if mod(i,2)==0 
        calcDist(j,1)=distances{i}(1,1);
        calcDist(j,2)=(sqrt(((distances{i}(1,2))-(distances{i}(2,2)))^2+((distances{i}(1,3))-(distances{i}(2,3)))^2))*pxl;
        j=j+1;
        else
        end
    end
    
    areas=[];
    j=1;
    for i=1:length(distances)
        if mod(i,2)==0 
        areas(j,1)=distances{i}(1,1);
        areas(j,2)=distances{i}(1,4);
        areas(j,3)=distances{i}(2,4);
        j=j+1;
        else
        end
    end
    
    perimeters=[];
    j=1;
    for i=1:length(distances)
        if mod(i,2)==0
        perimeters(j,1)=distances{i}(1,1);
        perimeters(j,2)=distances{i}(1,5);
        perimeters(j,3)=distances{i}(2,5);
        j=j+1;
        else
        end
    end
    
    circularity=[];
    j=1;
    for i=1:length(distances)
        if mod(i,2)==0
        circularity(j,1)=distances{i}(1,1);
        circularity(j,2)=distances{i}(1,6);
        circularity(j,3)=distances{i}(2,6);
        j=j+1;
        else
        end
    end
    
    deltaDist=[];
    directionalDist=[];
    for i=minFr:length(calcDist)
        if calcDist(i,1)<=maxFr+1
            if calcDist((i-1),1)==calcDist(i,1)-1
            deltaDist(i-1,1)=calcDist(i,1);
            deltaDist(i-1,2)=abs(calcDist(i,2)-calcDist((i-1),2));
            directionalDist(i-1,1)=calcDist(i,1);
            directionalDist(i-1,2)=(calcDist(i,2)-calcDist((i-1),2));
            else
            end
        else
        end
    end
    
    if isempty(deltaDist)~=1
        NZdeltaDist=nonzeros(deltaDist);
        NZdeltaDist=reshape(NZdeltaDist,size(NZdeltaDist,1)/size(deltaDist,2),size(deltaDist,2));
        NZdirectionalDist=nonzeros(directionalDist);
        NZdirectionalDist=reshape(NZdirectionalDist,size(NZdirectionalDist,1)/size(directionalDist,2),size(directionalDist,2));
        
        mincheck=[mincheck;min(NZdeltaDist(:,1))];
            
        allDist{k,5}(:,1) = NZdeltaDist(:,1);
        allDist{k,5}(:,2) = NZdeltaDist(:,2);
        allDist{k,6}(:,1) = NZdirectionalDist(:,1);
        allDist{k,6}(:,2) = NZdirectionalDist(:,2);
    else
    end
    

%     corrMotion=[];
%     for i=1:length(distances)
%         corrMotion(i,1)=distances{i}(1,1);
%         corrMotion(i,2)=
    
allDist{k,1}(:,1) = calcDist(:,1);
allDist{k,1}(:,2) = calcDist(:,2);
allDist{k,2}(:,1) = areas(:,1);
allDist{k,2}(:,2) = areas(:,2);
allDist{k,2}(:,3) = areas(:,3);
allDist{k,3}(:,1) = perimeters(:,1);
allDist{k,3}(:,2) = perimeters(:,2);
allDist{k,3}(:,3) = perimeters(:,3);
allDist{k,4}(:,1) = circularity(:,1);
allDist{k,4}(:,2) = circularity(:,2);
allDist{k,4}(:,3) = circularity(:,3);

allDist{k,7}(1,:) = splitName{1,1};


    
end


save(analysisTitle)

close all

min(mincheck)


iii=1;

% for i=1:size(allDist,1)
% %     max(allDist{i,5}(:,2))
%     for j=1:length(allDist{i,5})
%         if allDist{i,5}(j,2)>=5
%             checkDIST{iii,1}=allDist{i,7};
%             checkDIST{iii,2}=allDist{i,6}(j,1);
%             checkDIST{iii,3}=allDist{i,6}(j,2);
%             iii=iii+1;
%         else
%         end
%     end
% end

toc

%%USEFUL PLOTTING CODE%%

% for i=1:length(allDist)
% hold on
% scatter(allDist{i,1}(:,1), allDist{i,1}(:,2),'.')
% end

% k=figure;
% hold on
% scatter(calcDist(:,1), calcDist(:,2),'.')

% for i=1:length(allDist)
% hold on
% scatter(allDist{i,1}(:,1), allDist{i,1}(:,2),'.')
% end

%%DISTANCES
% figure
% hold on
% scatter(allDist{1,1}(:,1), allDist{1,1}(:,2),'.','bl')
% scatter(allDist{2,1}(:,1), allDist{2,1}(:,2),'.','m')
% scatter(allDist{3,1}(:,1), allDist{3,1}(:,2),'.','m')
% scatter(allDist{4,1}(:,1), allDist{4,1}(:,2),'.','m')
% scatter(allDist{5,1}(:,1), allDist{5,1}(:,2),'.','m')
% scatter(allDist{6,1}(:,1), allDist{6,1}(:,2),'.','m')
% scatter(allDist{8,1}(:,1), allDist{8,1}(:,2),'.','m')
% scatter(allDist{9,1}(:,1), allDist{9,1}(:,2),'.','m')
% scatter(allDist{10,1}(:,1), allDist{10,1}(:,2),'.','m')
% scatter(allDist{11,1}(:,1), allDist{11,1}(:,2),'.','m')
% scatter(allDist{12,1}(:,1), allDist{12,1}(:,2),'.','m')
% scatter(allDist{13,1}(:,1), allDist{13,1}(:,2),'.','m')
% scatter(allDist{14,1}(:,1), allDist{14,1}(:,2),'.','m')

%%AREAS
% figure
% hold on
% scatter(allDist{1,2}(:,1), allDist{1,2}(:,2),'.','bl')
% scatter(allDist{2,2}(:,1), allDist{2,2}(:,2),'.','m')
% scatter(allDist{3,2}(:,1), allDist{3,2}(:,2),'.','m')
% scatter(allDist{4,2}(:,1), allDist{4,2}(:,2),'.','m')
% scatter(allDist{5,2}(:,1), allDist{5,2}(:,2),'.','m')
% scatter(allDist{6,2}(:,1), allDist{6,2}(:,2),'.','m')
% scatter(allDist{7,2}(:,1), allDist{7,2}(:,2),'.','m')
% scatter(allDist{8,2}(:,1), allDist{8,2}(:,2),'.','m')
% scatter(allDist{9,2}(:,1), allDist{9,2}(:,2),'.','m')
% scatter(allDist{10,2}(:,1), allDist{10,2}(:,2),'.','m')
% scatter(allDist{11,2}(:,1), allDist{11,2}(:,2),'.','m')
% scatter(allDist{12,2}(:,1), allDist{12,2}(:,2),'.','m')
% scatter(allDist{13,2}(:,1), allDist{13,2}(:,2),'.','m')
% scatter(allDist{14,2}(:,1), allDist{14,2}(:,2),'.','m')

