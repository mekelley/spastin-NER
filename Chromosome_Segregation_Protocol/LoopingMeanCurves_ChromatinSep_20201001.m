%clearvars -except FifteenMinHeLa FifteenMinRPE1 FifteenMinWTDose FifteenMinN386CDose FifteenMinshSPAST
clear all
close all

% Specify the folder where the files live.
theFolder = 'D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_Chromosome_Separation\CellLines_wekaAnalysis\MatLabAnalyses\HeLaP';
testTitle = '20210302_HeLaP';
cd('D:\Rockefeller\PhD_Research_Kapoor\Publication\MEK_LC_2020\Data\Live_Chromosome_Separation\CellLines_wekaAnalysis\MatLabAnalyses\Analysis');

d1spz=['1 ',char(181),'M spastazoline'];
d2p5spz=['2.5 ',char(181),'M spastazoline'];
d5spz=['5 ',char(181),'M spastazoline'];
d10spz=['10 ',char(181),'M spastazoline'];
wtspz=strcat('WT',{' '},d1spz);
N386Cspz=strcat('N386C',{' '},d1spz);
micron=[char(181),'m'];

% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(theFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', theFolder);
  uiwait(warndlg(errorMessage));
  return;
end

% Get a list of all files in the folder with the desired file name
% pattern.c
fPattern = fullfile(theFolder, '*.mat'); % Change to whatever pattern you need.
matFiles = dir(fPattern);

FifteenMin={};
deltaD={};
directionalD={};

%for n = 5:8
for n = 1:length(matFiles) 
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp

    %baseName = matFiles((n-4)).name;
    baseName = matFiles(n).name;
    fullName = fullfile(theFolder, baseName);
    [filepath,name,ext] = fileparts(fullName);
    matName=name;
    %spName = split(name,'_');
    
    load(fullfile(filepath, strcat(name,ext)));

    k=1;
    for i=1:size(allDist,1)
        AllData{k,1}=allDist{i,1};
        k=k+1;
    end

        BINsize=1; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,1));
          ShortestTime(t)=min(AllData{t}(:,1));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,1)>=TimeStamp(i) & AllData{ii}(:,1)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,2))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=nanmean(BinIntInt{m});
              stdIntInt(m)=nanstd(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,1}=matName;
    totalExps{n,2}=TimeStamp;
    totalExps{n,3}=meanIntInt;
    totalExps{n,4}=stdIntInt;
    totalExps{n,5}=SEMIntInt;
    
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp

    k=1;
    for i=1:size(allDist,1)
        AllData{k,1}=allDist{k,2}(:,1:2:3);
        k=k+1;
    end

        BINsize=1; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,1));
          ShortestTime(t)=min(AllData{t}(:,1));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,1)>=TimeStamp(i) & AllData{ii}(:,1)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,2))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=nanmean(BinIntInt{m});
              stdIntInt(m)=nanstd(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,6}=TimeStamp;
    totalExps{n,7}=meanIntInt;
    totalExps{n,8}=stdIntInt;
    totalExps{n,9}=SEMIntInt;
    
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp

    k=1;
    for i=1:size(allDist,1)
        AllData{k,1}=allDist{k,2}(:,1:1:2);
        k=k+1;
    end

        BINsize=1; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,1));
          ShortestTime(t)=min(AllData{t}(:,1));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,1)>=TimeStamp(i) & AllData{ii}(:,1)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,2))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=nanmean(BinIntInt{m});
              stdIntInt(m)=nanstd(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,10}=TimeStamp;
    totalExps{n,11}=meanIntInt;
    totalExps{n,12}=stdIntInt;
    totalExps{n,13}=SEMIntInt;
    
    %%added%%
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp

    k=1;
    for i=1:size(allDist,1)
        AllData{k,1}=allDist{k,4}(:,1:2:3);
        k=k+1;
    end

        BINsize=1; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,1));
          ShortestTime(t)=min(AllData{t}(:,1));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,1)>=TimeStamp(i) & AllData{ii}(:,1)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,2))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=nanmean(BinIntInt{m});
              stdIntInt(m)=nanstd(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,14}=TimeStamp;
    totalExps{n,15}=meanIntInt;
    totalExps{n,16}=stdIntInt;
    totalExps{n,17}=SEMIntInt;
    
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp

    k=1;
    for i=1:size(allDist,1)
        AllData{k,1}=allDist{k,4}(:,1:1:2);
        k=k+1;
    end

        BINsize=1; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,1));
          ShortestTime(t)=min(AllData{t}(:,1));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,1)>=TimeStamp(i) & AllData{ii}(:,1)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,2))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=nanmean(BinIntInt{m});
              stdIntInt(m)=nanstd(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,18}=TimeStamp;
    totalExps{n,19}=meanIntInt;
    totalExps{n,20}=stdIntInt;
    totalExps{n,21}=SEMIntInt;
    
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp
    
    k=1;
    for i=1:size(allDist,1)
        AllData{k,1}=allDist{i,5};
        k=k+1;
    end

        BINsize=1; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,1));
          ShortestTime(t)=min(AllData{t}(:,1));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,1)>=TimeStamp(i) & AllData{ii}(:,1)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,2))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=nanmean(BinIntInt{m});
              stdIntInt(m)=nanstd(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,22}=TimeStamp;
    totalExps{n,23}=meanIntInt;
    totalExps{n,24}=stdIntInt;
    totalExps{n,25}=SEMIntInt;
    
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp

        k=1;
    for i=1:size(allDist,1)
        AllData{k,1}=allDist{i,6};
        k=k+1;
    end

        BINsize=1; %interval size (s)

    for t=1:length(AllData)
          LongestTime(t)=max(AllData{t}(:,1));
          ShortestTime(t)=min(AllData{t}(:,1));
      end
          MaxTimeStamp=max(LongestTime);
          MinTimeStamp=min(ShortestTime);

          TimeStamp=MinTimeStamp:BINsize:MaxTimeStamp;

          %make bins for integrated intensity values
          BinIntInt={};
             for b=1:length(TimeStamp)
                  BinIntInt{b}=[];
             end

          for i=1:length(BinIntInt)-1
              for ii=1:length(AllData)
                  idx=find(AllData{ii}(:,1)>=TimeStamp(i) & AllData{ii}(:,1)<=TimeStamp(i+1));
                  if isempty(idx)==0
                      BinIntInt{i}=[BinIntInt{i}; (AllData{ii}(idx,2))];
                  else                
                  clear idx
                end
              end
          end
          for m=1:length(BinIntInt)
              meanIntInt(m)=nanmean(BinIntInt{m});
              stdIntInt(m)=nanstd(BinIntInt{m});
              SEMIntInt(m)=stdIntInt(m)./length(BinIntInt{m});
          end
          
    totalExps{n,26}=TimeStamp;
    totalExps{n,27}=meanIntInt;
    totalExps{n,28}=stdIntInt;
    totalExps{n,29}=SEMIntInt;
    
    clear AllData meanIntInt stdIntInt SEMIntInt BinIntInt LongestTime ShortestTime MaxTimeStamp MinTimeStamp TimeStamp
    
for i=1:size(allDist,1)
    for j=1:length(allDist{i,1})
        if allDist{i,1}(j)==15
            FifteenMin{i,n}=allDist{i,1}(j,2);
        else
        end
    end
end


for i=1:size(allDist,1)
    for j=1:length(allDist{i,1})
        if allDist{i,1}(j)==12
            TwelveMin{i,n}=allDist{i,1}(j,2);
        else
        end
    end
end


for i=1:size(allDist,1)
    for j=1:length(allDist{i,1})
        if allDist{i,1}(j)==10
            TenMin{i,n}=allDist{i,1}(j,2);
        else
        end
    end
end

for i=1:size(allDist,1)
    for j=1:length(allDist{i,1})
        if allDist{i,1}(j)==8
            EightMin{i,n}=allDist{i,1}(j,2);
        else
        end
    end
end

for i=1:size(allDist,1)
    for j=1:length(allDist{i,1})
        if allDist{i,1}(j)==5
            FiveMin{i,n}=allDist{i,1}(j,2);
        else
        end
    end
end

for i=1:size(allDist,1)
    for j=1:length(allDist{i,1})
        if allDist{i,1}(j)==3
            ThreeMin{i,n}=allDist{i,1}(j,2);
        else
        end
    end
end



for i=1:size(allDist,1)
        deltaD{i,n}=allDist{i,5}(:,2);
end
    
for i=1:size(allDist,1)
        directionalD{i,n}=allDist{i,6}(:,2);
end
    
end

save(testTitle);


%%%%%%%%%%%%%%%%%%%%%2 items %%%%%%%%%%%%%%%%%%%%
deltaDMSO=[];
deltaSPZ=[];
szdelta=size(deltaD);
for i=1:szdelta(2)
    for j=1:szdelta(1)
        if i==1
        deltaDMSO=[deltaDMSO;deltaD{j,i}];
        elseif i==2
            deltaSPZ=[deltaSPZ;deltaD{j,i}];
        end
    end
end

% %%%%%%%%%%%%%%%%%%%%4 items dose%%%%%%%%%%%%
% deltaDMSO=[];
% delta2p5=[];
% delta5=[];
% delta10=[];
% szdelta=size(deltaD);
% for i=1:szdelta(2)
%     for j=1:szdelta(1)
%         if i==1
%         deltaDMSO=[deltaDMSO;deltaD{j,i}];
%         elseif i==2
%             delta2p5=[delta2p5;deltaD{j,i}];
%         elseif i==3
%             delta5=[delta5;deltaD{j,i}];
%         elseif i==4
%             delta10=[delta10;deltaD{j,i}];
%         end
%     end
% end

% %%%%%%%%%%%%%%%%%%%%4 items shSPAST%%%%%%%%%%%%
% deltaCtrlDMSO=[];
% deltaCtrlSPZ=[];
% deltaDoxDMSO=[];
% deltaDoxSPZ=[];
% szdelta=size(deltaD);
% for i=1:szdelta(2)
%     for j=1:szdelta(1)
%         if i==1
%         deltaCtrlDMSO=[deltaCtrlDMSO;deltaD{j,i}];
%         elseif i==2
%             deltaCtrlSPZ=[deltaCtrlSPZ;deltaD{j,i}];
%         elseif i==3
%             deltaDoxDMSO=[deltaDoxDMSO;deltaD{j,i}];
%         elseif i==4
%             deltaDoxSPZ=[deltaDoxSPZ;deltaD{j,i}];
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%
%  
% % %PLOTS: experiments must be in DMSO then SPZ order for each celltype in totalExps 
% % %for DMSO to be in black and spastazoline in magenta
% % for n=1:length(matFiles)
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         figTitle=split(totalExps{n,1}, '_');
% %         %plot centroid distances
% %         figure;
% %         subplot(1,3,1)
% %         hold on
% % %         title('Centroid Distance')
% %         shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'bl', 1);
% %         shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'm', 1);
% %         xlabel('Time (min)') 
% %         ylabel('Distance (microns)')
% %         axis([0 100 0 70])
% %         %plot areas of nuclei I
% %         subplot(1,3,2)
% %         hold on
% %         %title('Nuclei Area I')
% %         shadedErrorBar(totalExps{n,6},totalExps{n,7},totalExps{n,8}, 'bl', 1); 
% %         shadedErrorBar(totalExps{x,6},totalExps{x,7},totalExps{x,8}, 'm', 1);
% %         xlabel('Time (min)') 
% %         ylabel('Area I (microns^2)')
% %         axis([0 100 0 250])
% %         %plot areas of nuclei II
% %         subplot(1,3,3)
% %         hold on
% %         %title('Nuclei Area II')
% %         shadedErrorBar(totalExps{n,10},totalExps{n,11},totalExps{n,12}, 'bl', 1);
% %         shadedErrorBar(totalExps{x,10},totalExps{x,11},totalExps{x,12}, 'm', 1);
% %         xlabel('Time (min)') 
% %         ylabel('Area II (microns^2)')
% %         mtit(strcat(char(figTitle(2)), ' Cells'))
% %     else
% %     end
% % end
% 
% % %PLOTS: experiments must be in DMSO then SPZ order for each celltype in totalExps 
% % %for DMSO to be in black and spastazoline in magenta
% % 
% % xmin=4;
% % ylim=70;
% % xlim=90;
% % yinset=20;
% % xinset=25;
% % ymin=8;
% % 
% % for n=1:length(matFiles)
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         figTitle=split(totalExps{n,1}, '_');
% %         %plot centroid distances
% %         h1=figure;
% % %         subplot(1,3,1)
% %         hold on
% % %         title('Centroid Distance')
% %         p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
% %         p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'b', 1);
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel('Internuclei distance (microns)')
% %         axis([xmin xlim 0 ylim])
% %         title(strcat(char(figTitle(2)), ' Cells'))
% %         legend([p1.patch,p2.patch],'0.1% DMSO', d10spz,'Location','northeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% %         set(h1, 'Position',[105   647   700   500]);
% %         MagInset(h1, -1, [xmin xinset ymin yinset], [10 50 40 70], {'NW','NW';'SE','SE'});
% % %       gp = plotboxpos(gca);
% % %       axes('position',[gp(1)+(.08)*gp(3) gp(2)+(.52)*gp(4) .38*gp(3) .48*gp(4)])
% %         %axes('Position',[.18 .55 .3 .5])
% % %         box on
% % %         hold on
% % %         shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
% % %         shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'b', 1);
% % %         axis([xmin xinset 0 yinset])
% %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% % end
% 
% %%%%%%DOSE RESPONSE%%%%%%%%%% 
% %PLOTS: experiments must be in DMSO then SPZ order for each celltype in totalExps 
% %for DMSO to be in black and spastazoline in magenta
% % 
% % xmin=4;
% % ylim=21;
% % xlim=20;
% % yinset=20;
% % xinset=20;
% % ymin=8;
% % 
% % n=5;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         figTitle=('HeLa N386C Cells');
% %         h1=figure;
% %         hold on
% %         p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
% %         p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'g', 1);
% %         p3=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'm', 1);
% %         p4=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'b', 1);
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel('Internuclei distance (microns)')
% %         axis([xmin xlim 0 ylim])
% %         title('HeLa N386C Cells')
% %         legend([p1.patch,p2.patch,p3.patch,p4.patch],'0.1% DMSO', d2p5spz, d5spz, d10spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% %     
% % %%%%%%DMSO vs 10uM SPZ%%%%%%%%% 
% % 
% % xmin=4;
% % ylim=21;
% % xlim=20;
% % yinset=20;
% % xinset=20;
% % ymin=8;
% % 
% % n=1;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         nn=n+4;
% %         xx=x+4;
% %         yy=y+4;
% %         zz=z+4;
% %         figTitle=strcat('Normalized DMSO');
% %         h1=figure;
% %         hold on
% %         p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
% %         p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'b', 1);
% % %         p3=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'm', 1);
% % %         p4=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'b', 1);
% %  %        p5=shadedErrorBar(totalExps{nn,2},totalExps{nn,3},totalExps{nn,4}, 'm', 1);
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel('Internuclei distance (microns)')
% %         axis([xmin xlim 0 ylim])
% %        % title('0.1% DMSO')
% %         legend([p1.patch,p2.patch],'0.1% DMSO', d10spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% %     
% % %%%%%%%%%%%%%NORMALIZATION%%%%%%%%%%%%%%%
% % for kk=1:4
% %     for mm=4:60
% %     	jj=find(totalExps{1,2}==mm);
% %         ll=find(totalExps{kk,2}==mm);
% %         normExps{kk,1}=totalExps{kk,1};
% %         normExps{kk,2}(mm-3)= mm;
% %         normExps{kk,3}(mm-3)=totalExps{kk,3}(ll)/totalExps{1,3}(jj);
% %     end
% % end
% % for kk=5:8
% %     for mm=4:60
% %     	jj=find(totalExps{5,2}==mm);
% %         ll=find(totalExps{kk,2}==mm);
% %         normExps{kk,1}=totalExps{kk,1};
% %         normExps{kk,2}(mm-3) = mm;
% %         normExps{kk,3}(mm-3)=totalExps{kk,3}(ll)/totalExps{5,3}(jj);
% %     end
% % end
% % 
% % 
% % xmin=4;
% % ylim=1.05;
% % xlim=20;
% % yinset=20;
% % xinset=20;
% % ymin=8;
% % 
% % n=1;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         figTitle=('HeLa WT Cells');
% %         h1=figure;
% %         hold on
% %         plot(normExps{n,2},normExps{n,3},'--k');
% %         plot(normExps{x,2},normExps{x,3},'-g');
% %         plot(normExps{y,2},normExps{y,3},'-m');
% %         plot(normExps{z,2},normExps{z,3},'-b');
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel({'Internuclei distance';'(fraction of control)'})
% %         axis([xmin xlim 0 ylim])
% %         title('HeLa WT Cells')
% %         legend('0.1% DMSO', d2p5spz, d5spz, d10spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% % 
% %     n=5;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         figTitle=('HeLa N386C Cells');
% %         h1=figure;
% %         hold on
% %         plot(normExps{n,2},normExps{n,3},'--k');
% %         plot(normExps{x,2},normExps{x,3},'-g');
% %         plot(normExps{y,2},normExps{y,3},'-m');
% %         plot(normExps{z,2},normExps{z,3},'-b');
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel({'Internuclei distance';'(fraction of control)'})
% %         axis([xmin xlim 0 ylim])
% %         title('HeLa N386C Cells')
% %         legend('0.1% DMSO', d2p5spz, d5spz, d10spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% %     
% %     %%%%%%%%%%shSPAST PLOTS%%%%%%%%%%%%%%%%%%%%%
% %     
% %     xmin=4;
% % ylim=21;
% % xlim=20;
% % yinset=20;
% % xinset=20;
% % ymin=8;
% % 
% % n=1;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         nn=n+4;
% %         xx=x+4;
% %         yy=y+4;
% %         zz=z+4;
% %         figTitle=strcat('shSPAST');
% %         h1=figure;
% %         hold on
% %       %p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
% %        % p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'm', 1);
% %        p3=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'b', 1);
% %         p4=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'g', 1);
% %         %p5=shadedErrorBar(totalExps{nn,2},totalExps{nn,3},totalExps{nn,4}, 'm', 1);
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel('Internuclei distance (microns)')
% %         axis([xmin xlim 0 ylim])
% % %         title('0.1% DMSO')
% %         %legend([p1.patch,p2.patch,p3.patch,p4.patch],'HeLa + 0.1% DMSO', '(-)shSPAST no dox', '(+)shSPAST + 0.1% DMSO', '(+)shSPAST + spastazoline','Location','southeast')
% %         legend([p3.patch,p4.patch],'HeLa shSPAST (+)','HeLa shSPAST (+) spastazoline','Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% %     
% % %%%%%%%%%%%%%%%%%%GET MEAN STD%%%%%%%%%%%%%%%
% % %%DMSO or 1 at time 15 min
% % totalExps{1,3}(13) %%mean
% % totalExps{1,4}(13) %%std
% % 
% % %%SPZ or 2 at time 15 min
% % totalExps{2,3}(13) %%mean
% % totalExps{2,4}(13) %%std
% % 
% % %%3 at time 15 min
% % totalExps{3,3}(13) %%mean
% % totalExps{3,4}(13) %%std
% % 
% % %%4 at time 15 min
% % totalExps{4,3}(13) %%mean
% % totalExps{4,4}(13) %%std
% % 
% % %%5 at time 15 min
% % totalExps{5,3}(13) %%mean
% % totalExps{5,4}(13) %%std
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %%%%%% 4h TX DMSO vs 1uM SPZ%%%%%%%%% 
% % 
% % xmin=4;
% % ylim=21;
% % xlim=20;
% % yinset=20;
% % xinset=20;
% % ymin=8;
% % 
% % n=1;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         nn=n+4;
% %         xx=x+4;
% %         yy=y+4;
% %         zz=z+4;
% %         figTitle=strcat('WT');
% %         h1=figure;
% %         hold on
% %         p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
% %         p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'b', 1);
% % %         p3=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'm', 1);
% % %         p4=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'b', 1);
% %  %        p5=shadedErrorBar(totalExps{nn,2},totalExps{nn,3},totalExps{nn,4}, 'm', 1);
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel('Internuclei distance (microns)')
% %         axis([xmin xlim 0 ylim])
% %         title('WT')
% %         legend([p1.patch,p2.patch],'0.1% DMSO', d1spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% %     
% %     n=1;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         nn=n+4;
% %         xx=x+4;
% %         yy=y+4;
% %         zz=z+4;
% %         figTitle=strcat('N386C');
% %         h1=figure;
% %         hold on
% %  %       p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
% %   %      p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'b', 1);
% %          p3=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'k', 1);
% %          p4=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'b', 1);
% %  %        p5=shadedErrorBar(totalExps{nn,2},totalExps{nn,3},totalExps{nn,4}, 'm', 1);
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel('Internuclei distance (microns)')
% %         axis([xmin xlim 0 ylim])
% %         title('N386C')
% %         legend([p3.patch,p4.patch],'0.1% DMSO', d1spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% %     
% %         n=1;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         nn=n+4;
% %         xx=x+4;
% %         yy=y+4;
% %         zz=z+4;
% %         h1=figure;
% %         hold on
% %         p1=plot(totalExps{n,2},totalExps{n,3}, '--k');
% %         p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'b', 1);
% %   %       p3=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'k', 1);
% %          p4=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'm', 1);
% %  %        p5=shadedErrorBar(totalExps{nn,2},totalExps{nn,3},totalExps{nn,4}, 'm', 1);
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel('Internuclei distance (microns)')
% %         axis([xmin xlim 0 ylim])
% %         title('Combined')
% %         legend([p1,p2.patch,p4.patch],['0.1% DMSO', wtspz, N386Cspz],'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% %     
% %     
% % %%%%%%%%%%%%%NORMALIZATION 4h TX%%%%%%%%%%%%%%%
% % for kk=1:4
% %     for mm=4:60
% %     	jj=find(totalExps{1,2}==mm);
% %         ll=find(totalExps{kk,2}==mm);
% %         normExps{kk,1}=totalExps{kk,1};
% %         normExps{kk,2}(mm-3)= mm;
% %         normExps{kk,3}(mm-3)=totalExps{kk,3}(ll)/totalExps{1,3}(jj);
% %     end
% % end
% % for kk=5:8
% %     for mm=4:60
% %     	jj=find(totalExps{5,2}==mm);
% %         ll=find(totalExps{kk,2}==mm);
% %         normExps{kk,1}=totalExps{kk,1};
% %         normExps{kk,2}(mm-3) = mm;
% %         normExps{kk,3}(mm-3)=totalExps{kk,3}(ll)/totalExps{5,3}(jj);
% %     end
% % end
% % 
% % 
% % xmin=4;
% % ylim=1.05;
% % xlim=20;
% % yinset=20;
% % xinset=20;
% % ymin=8;
% % 
% % n=1;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         figTitle=('HeLa WT Cells');
% %         h1=figure;
% %         hold on
% %         plot(normExps{n,2},normExps{n,3},'--k');
% %         plot(normExps{x,2},normExps{x,3},'-g');
% %         plot(normExps{y,2},normExps{y,3},'-m');
% %         plot(normExps{z,2},normExps{z,3},'-b');
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel({'Internuclei distance';'(fraction of control)'})
% %         axis([xmin xlim 0 ylim])
% %         title('HeLa WT Cells')
% %         legend('0.1% DMSO', d2p5spz, d5spz, d10spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
% % 
% %     n=5;
% %     if mod(n,2)~=0 
% %         x=n+1;
% %         y=n+2;
% %         z=n+3;
% %         figTitle=('HeLa N386C Cells');
% %         h1=figure;
% %         hold on
% %         plot(normExps{n,2},normExps{n,3},'--k');
% %         plot(normExps{x,2},normExps{x,3},'-g');
% %         plot(normExps{y,2},normExps{y,3},'-m');
% %         plot(normExps{z,2},normExps{z,3},'-b');
% %         xlabel('Time post-anaphase onset (min)') 
% %         ylabel({'Internuclei distance';'(fraction of control)'})
% %         axis([xmin xlim 0 ylim])
% %         title('HeLa N386C Cells')
% %         legend('0.1% DMSO', d2p5spz, d5spz, d10spz,'Location','southeast')
% %         legend('boxoff')
% %         set(gca,'FontSize', 16);
% % %         set(h1, 'Position',[105   647   700   500]);
% % %         MagInset(h1, -1, [6 xinset ymin yinset], [8 35 40 70], {'NW','NW';'SE','SE'});
% % %         set(gca,'FontSize', 16);
% %         
% %     else
% %     end
%     
    %%%%%%%More plotting%%%%%%%%%%%%%
figure('Name',testTitle,'units','normalized','outerposition',[0 0 1 1]);
hold on;

xmin=min(mincheck);
ylim=40;
xlim=maxFr;
ymin=0;

n=1;
    if mod(n,2)~=0 
        x=n+1;
        subplot(2,4,[1 5]);
        figTitle=strcat('Centroids');
        title(figTitle);
        %h1=figure;
        hold on
        p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
        p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'm', 1);
        xlabel('time post-anaphase onset (min)') 
        ylabel(strcat('internuclei distance (',micron,')'))
        axis([xmin xlim 0 ylim])
        legend([p1.patch,p2.patch],'0.1% DMSO',d10spz,'Location','southeast')
        legend('boxoff')
        set(gca,'FontSize', 16);

        
    else
    end


    
% xmin=5;
ylim=200;
% xlim=60;
ymin=0;

    n=1;
    if mod(n,2)~=0 
        x=n+1;
        figTitle=strcat('Area I');
        %h2=figure;
        subplot(2,4,2);
        title(figTitle);
        hold on
        p3=shadedErrorBar(totalExps{n,2},totalExps{n,7},totalExps{n,8}, 'k', 1);
        p4=shadedErrorBar(totalExps{x,2},totalExps{x,7},totalExps{x,8}, 'm', 1);
        xlabel('time post-anaphase onset (min)') 
        ylabel(strcat('nuclei area (',micron,'^{2})'))
        axis([xmin xlim 0 ylim])
        legend([p3.patch,p4.patch],'0.1% DMSO',d10spz,'Location','southeast')
        legend('boxoff')
        set(gca,'FontSize', 16);

        
    else
    end
        
    n=1;
    if mod(n,2)~=0 
        x=n+1;
        figTitle=strcat('Area II');
        %h3=figure;
        subplot(2,4,6);
        title(figTitle);
        hold on
        p5=shadedErrorBar(totalExps{n,2},totalExps{n,11},totalExps{n,12}, 'k', 1);
        p6=shadedErrorBar(totalExps{x,2},totalExps{x,11},totalExps{x,12}, 'm', 1);
        xlabel('time post-anaphase onset (min)') 
        ylabel(strcat('nuclei area (',micron,'^{2})'))
        axis([xmin xlim 0 ylim])
        legend([p5.patch,p6.patch],'0.1% DMSO',d10spz,'Location','southeast')
        legend('boxoff')
        set(gca,'FontSize', 16);

        
    else
    end
    
% xmin=5;
ylim=1;
% xlim=60;
ymin=0;

    n=1;
    if mod(n,2)~=0 
        x=n+1;
        figTitle=strcat('Circularity I');
        %h4=figure;
        subplot(2,4,3);
        title(figTitle);
        hold on
        p7=shadedErrorBar(totalExps{n,2},totalExps{n,15},totalExps{n,16}, 'k', 1);
        p8=shadedErrorBar(totalExps{x,2},totalExps{x,15},totalExps{x,16}, 'm', 1);
        xlabel('time post-anaphase onset (min)') 
        ylabel(strcat('nuclei circularity'))
        axis([xmin xlim 0 ylim])
        legend([p7.patch,p8.patch],'0.1% DMSO',d10spz,'Location','southeast')
        legend('boxoff')
        set(gca,'FontSize', 16);

        
    else
    end
        
    n=1;
    if mod(n,2)~=0 
        x=n+1;
        figTitle=strcat('Circularity II');
        %h5=figure;
        subplot(2,4,7);
        title(figTitle);
        hold on
        p9=shadedErrorBar(totalExps{n,2},totalExps{n,19},totalExps{n,20}, 'k', 1);
        p10=shadedErrorBar(totalExps{x,2},totalExps{x,19},totalExps{x,20}, 'm', 1);
        xlabel('time post-anaphase onset (min)') 
        ylabel(strcat('nuclei circularity'))
        axis([xmin xlim 0 ylim])
        legend([p9.patch,p10.patch],'0.1% DMSO',d10spz,'Location','southeast')
        legend('boxoff')
        set(gca,'FontSize', 16);

        
    else
    end
    
% xmin=5;
ylim=1.5;
% xlim=60;
ymin=0;

    n=1;
    if mod(n,2)~=0 
        x=n+1;
        figTitle=strcat('Delta Distance');
        %h6=figure;
        subplot(2,4,4);
        title(figTitle);
        hold on
        p11=shadedErrorBar(totalExps{n,22},totalExps{n,23},totalExps{n,24}, 'k', 1);
        p12=shadedErrorBar(totalExps{x,22},totalExps{x,23},totalExps{x,24}, 'm', 1);
        xlabel('time post-anaphase onset (min)') 
        ylabel(strcat('\Delta distance between centroids (',micron,')'))
        axis([xmin xlim 0 ylim])
        legend([p11.patch,p12.patch],'0.1% DMSO',d10spz,'Location','northwest')
        legend('boxoff')
        set(gca,'FontSize', 16);

        
    else
    end
 

% xmin=5;
ylim=1.5;
% xlim=60;
ymin=-1;

    n=1;
    if mod(n,2)~=0 
        x=n+1;
        figTitle=strcat('Delta Distance Direction');
        %h7=figure;
        subplot(2,4,8);
        title(figTitle);
        hold on
        p13=shadedErrorBar(totalExps{n,26},totalExps{n,27},totalExps{n,28}, 'k', 1);
        p14=shadedErrorBar(totalExps{x,26},totalExps{x,27},totalExps{x,28}, 'm', 1);
        plot([xmin,xlim],[0,0],'k--')
        xlabel('time post-anaphase onset (min)') 
        ylabel(strcat('\Delta distance between centroids (',micron,')'))
        axis([xmin xlim ymin ylim])
        legend([p13.patch,p14.patch],'0.1% DMSO',d10spz,'Location','northwest')
        legend('boxoff')
        set(gca,'FontSize', 16);

        
    else
    end
    
savefig(testTitle);

%     %%%%%%%More dose plotting%%%%%%%%%%%%%
% figure('Name',testTitle,'units','normalized','outerposition',[0 0 1 1]);
% hold on;
% 
% xmin=5;
% ylim=40;
% xlim=60;
% ymin=0;
% 
% n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         subplot(2,4,[1 5]);
%         figTitle=strcat('Centroids');
%         title(figTitle);
%         %h1=figure;
%         hold on
%         p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
%         p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'g', 1);
%         p1a=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'm', 1);
%         p2b=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('internuclei distance (',micron,')'))
%         axis([xmin xlim 0 ylim])
%         legend([p1.patch,p2.patch,p1a.patch,p2b.patch],'0.1% DMSO',d2p5spz,d5spz,d10spz,'Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
% 
% 
%     
% xmin=5;
% ylim=200;
% xlim=60;
% ymin=0;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Area I');
%         %h2=figure;
%         subplot(2,4,2);
%         title(figTitle);
%         hold on
%         p3=shadedErrorBar(totalExps{n,2},totalExps{n,7},totalExps{n,8}, 'k', 1);
%         p4=shadedErrorBar(totalExps{x,2},totalExps{x,7},totalExps{x,8}, 'g', 1);
%         p3a=shadedErrorBar(totalExps{y,2},totalExps{y,7},totalExps{y,8}, 'm', 1);
%         p4a=shadedErrorBar(totalExps{z,2},totalExps{z,7},totalExps{z,8}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei area (',micron,'^{2})'))
%         axis([xmin xlim 0 ylim])
%         legend([p3.patch,p4.patch,p3a.patch,p4a.patch],'0.1% DMSO',d2p5spz,d5spz,d10spz,'Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%         
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Area II');
%         %h3=figure;
%         subplot(2,4,6);
%         title(figTitle);
%         hold on
%         p5=shadedErrorBar(totalExps{n,2},totalExps{n,11},totalExps{n,12}, 'k', 1);
%         p6=shadedErrorBar(totalExps{x,2},totalExps{x,11},totalExps{x,12}, 'g', 1);
%         p5a=shadedErrorBar(totalExps{y,2},totalExps{y,11},totalExps{y,12}, 'm', 1);
%         p6a=shadedErrorBar(totalExps{z,2},totalExps{z,11},totalExps{z,12}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei area (',micron,'^{2})'))
%         axis([xmin xlim 0 ylim])
%         legend([p5.patch,p6.patch,p5a.patch,p6a.patch],'0.1% DMSO',d2p5spz,d5spz,d10spz,'Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% xmin=5;
% ylim=1;
% xlim=60;
% ymin=0;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Circularity I');
%         %h4=figure;
%         subplot(2,4,3);
%         title(figTitle);
%         hold on
%         p7=shadedErrorBar(totalExps{n,2},totalExps{n,15},totalExps{n,16}, 'k', 1);
%         p8=shadedErrorBar(totalExps{x,2},totalExps{x,15},totalExps{x,16}, 'g', 1);
%         p7a=shadedErrorBar(totalExps{y,2},totalExps{y,15},totalExps{y,16}, 'm', 1);
%         p8a=shadedErrorBar(totalExps{z,2},totalExps{z,15},totalExps{z,16}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei circularity'))
%         axis([xmin xlim 0 ylim])
%         legend([p7.patch,p8.patch,p7a.patch,p8a.patch],'0.1% DMSO',d2p5spz,d5spz,d10spz,'Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%         
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Circularity II');
%         %h5=figure;
%         subplot(2,4,7);
%         title(figTitle);
%         hold on
%         p9=shadedErrorBar(totalExps{n,2},totalExps{n,19},totalExps{n,20}, 'k', 1);
%         p10=shadedErrorBar(totalExps{x,2},totalExps{x,19},totalExps{x,20}, 'g', 1);
%         p9a=shadedErrorBar(totalExps{y,2},totalExps{y,19},totalExps{y,20}, 'm', 1);
%         p10a=shadedErrorBar(totalExps{z,2},totalExps{z,19},totalExps{z,20}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei circularity'))
%         axis([xmin xlim 0 ylim])
%         legend([p9.patch,p10.patch,p9a.patch,p10a.patch],'0.1% DMSO',d2p5spz,d5spz,d10spz,'Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% xmin=5;
% ylim=3;
% xlim=60;
% ymin=0;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Delta Distance');
%         %h6=figure;
%         subplot(2,4,4);
%         title(figTitle);
%         hold on
%         p11=shadedErrorBar(totalExps{n,22},totalExps{n,23},totalExps{n,24}, 'k', 1);
%         p12=shadedErrorBar(totalExps{x,22},totalExps{x,23},totalExps{x,24}, 'g', 1);
%         p11a=shadedErrorBar(totalExps{y,22},totalExps{y,23},totalExps{y,24}, 'm', 1);
%         p12a=shadedErrorBar(totalExps{z,22},totalExps{z,23},totalExps{z,24}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('\Delta distance between centroids (',micron,')'))
%         axis([xmin xlim 0 ylim])
%         legend([p11.patch,p12.patch,p11a.patch,p12a.patch],'0.1% DMSO',d2p5spz,d5spz,d10spz,'Location','northwest')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% xmin=5;
% ylim=2.5;
% xlim=60;
% ymin=-2.5;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Delta Distance Direction');
%         %h7=figure;
%         subplot(2,4,8);
%         title(figTitle);
%         hold on
%         p13=shadedErrorBar(totalExps{n,26},totalExps{n,27},totalExps{n,28}, 'k', 1);
%         p14=shadedErrorBar(totalExps{x,26},totalExps{x,27},totalExps{x,28}, 'g', 1);
%         p13a=shadedErrorBar(totalExps{y,26},totalExps{y,27},totalExps{y,28}, 'm', 1);
%         p14a=shadedErrorBar(totalExps{z,26},totalExps{z,27},totalExps{z,28}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('\Delta distance between centroids (',micron,')'))
%         axis([xmin xlim ymin ylim])
%         legend([p13.patch,p14.patch,p13a.patch,p14a.patch],'0.1% DMSO',d2p5spz,d5spz,d10spz,'Location','northwest')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% savefig(testTitle);

%     %%%%%%%More shSPAST plotting%%%%%%%%%%%%%
% figure('Name',testTitle,'units','normalized','outerposition',[0 0 1 1]);
% hold on;
% 
% xmin=5;
% ylim=40;
% xlim=60;
% ymin=0;
% 
% n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         subplot(2,4,[1 5]);
%         figTitle=strcat('Centroids');
%         title(figTitle);
%         %h1=figure;
%         hold on
%         p1=shadedErrorBar(totalExps{n,2},totalExps{n,3},totalExps{n,4}, 'k', 1);
%         p2=shadedErrorBar(totalExps{x,2},totalExps{x,3},totalExps{x,4}, 'g', 1);
%         p1a=shadedErrorBar(totalExps{y,2},totalExps{y,3},totalExps{y,4}, 'm', 1);
%         p2b=shadedErrorBar(totalExps{z,2},totalExps{z,3},totalExps{z,4}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('internuclei distance (',micron,')'))
%         axis([xmin xlim 0 ylim])
%         legend([p1.patch,p2.patch,p1a.patch,p2b.patch],'shSPAST ctrl DMSO','shSPAST ctrl SPZ','shSPAST dox DMSO', 'shSPAST dox SPZ','Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
% 
% 
%     
% xmin=5;
% ylim=200;
% xlim=60;
% ymin=0;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Area I');
%         %h2=figure;
%         subplot(2,4,2);
%         title(figTitle);
%         hold on
%         p3=shadedErrorBar(totalExps{n,2},totalExps{n,7},totalExps{n,8}, 'k', 1);
%         p4=shadedErrorBar(totalExps{x,2},totalExps{x,7},totalExps{x,8}, 'g', 1);
%         p3a=shadedErrorBar(totalExps{y,2},totalExps{y,7},totalExps{y,8}, 'm', 1);
%         p4a=shadedErrorBar(totalExps{z,2},totalExps{z,7},totalExps{z,8}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei area (',micron,'^{2})'))
%         axis([xmin xlim 0 ylim])
%         legend([p3.patch,p4.patch,p3a.patch,p4a.patch],'shSPAST ctrl DMSO','shSPAST ctrl SPZ','shSPAST dox DMSO', 'shSPAST dox SPZ','Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%         
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Area II');
%         %h3=figure;
%         subplot(2,4,6);
%         title(figTitle);
%         hold on
%         p5=shadedErrorBar(totalExps{n,2},totalExps{n,11},totalExps{n,12}, 'k', 1);
%         p6=shadedErrorBar(totalExps{x,2},totalExps{x,11},totalExps{x,12}, 'g', 1);
%         p5a=shadedErrorBar(totalExps{y,2},totalExps{y,11},totalExps{y,12}, 'm', 1);
%         p6a=shadedErrorBar(totalExps{z,2},totalExps{z,11},totalExps{z,12}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei area (',micron,'^{2})'))
%         axis([xmin xlim 0 ylim])
%         legend([p5.patch,p6.patch,p5a.patch,p6a.patch],'shSPAST ctrl DMSO','shSPAST ctrl SPZ','shSPAST dox DMSO', 'shSPAST dox SPZ','Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% xmin=5;
% ylim=1;
% xlim=60;
% ymin=0;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Circularity I');
%         %h4=figure;
%         subplot(2,4,3);
%         title(figTitle);
%         hold on
%         p7=shadedErrorBar(totalExps{n,2},totalExps{n,15},totalExps{n,16}, 'k', 1);
%         p8=shadedErrorBar(totalExps{x,2},totalExps{x,15},totalExps{x,16}, 'g', 1);
%         p7a=shadedErrorBar(totalExps{y,2},totalExps{y,15},totalExps{y,16}, 'm', 1);
%         p8a=shadedErrorBar(totalExps{z,2},totalExps{z,15},totalExps{z,16}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei circularity'))
%         axis([xmin xlim 0 ylim])
%         legend([p7.patch,p8.patch,p7a.patch,p8a.patch],'shSPAST ctrl DMSO','shSPAST ctrl SPZ','shSPAST dox DMSO', 'shSPAST dox SPZ','Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%         
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Circularity II');
%         %h5=figure;
%         subplot(2,4,7);
%         title(figTitle);
%         hold on
%         p9=shadedErrorBar(totalExps{n,2},totalExps{n,19},totalExps{n,20}, 'k', 1);
%         p10=shadedErrorBar(totalExps{x,2},totalExps{x,19},totalExps{x,20}, 'g', 1);
%         p9a=shadedErrorBar(totalExps{y,2},totalExps{y,19},totalExps{y,20}, 'm', 1);
%         p10a=shadedErrorBar(totalExps{z,2},totalExps{z,19},totalExps{z,20}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('nuclei circularity'))
%         axis([xmin xlim 0 ylim])
%         legend([p9.patch,p10.patch,p9a.patch,p10a.patch],'shSPAST ctrl DMSO','shSPAST ctrl SPZ','shSPAST dox DMSO', 'shSPAST dox SPZ','Location','southeast')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% xmin=5;
% ylim=3;
% xlim=60;
% ymin=0;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Delta Distance');
%         %h6=figure;
%         subplot(2,4,4);
%         title(figTitle);
%         hold on
%         p11=shadedErrorBar(totalExps{n,22},totalExps{n,23},totalExps{n,24}, 'k', 1);
%         p12=shadedErrorBar(totalExps{x,22},totalExps{x,23},totalExps{x,24}, 'g', 1);
%         p11a=shadedErrorBar(totalExps{y,22},totalExps{y,23},totalExps{y,24}, 'm', 1);
%         p12a=shadedErrorBar(totalExps{z,22},totalExps{z,23},totalExps{z,24}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('\Delta distance between centroids (',micron,')'))
%         axis([xmin xlim 0 ylim])
%         legend([p11.patch,p12.patch,p11a.patch,p12a.patch],'shSPAST ctrl DMSO','shSPAST ctrl SPZ','shSPAST dox DMSO', 'shSPAST dox SPZ','Location','northwest')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% xmin=5;
% ylim=2.5;
% xlim=60;
% ymin=-2.5;
% 
%     n=1;
%     if mod(n,2)~=0 
%         x=n+1;
%         y=n+2;
%         z=n+3;
%         figTitle=strcat('Delta Distance Direction');
%         %h7=figure;
%         subplot(2,4,8);
%         title(figTitle);
%         hold on
%         p13=shadedErrorBar(totalExps{n,26},totalExps{n,27},totalExps{n,28}, 'k', 1);
%         p14=shadedErrorBar(totalExps{x,26},totalExps{x,27},totalExps{x,28}, 'g', 1);
%         p13a=shadedErrorBar(totalExps{y,26},totalExps{y,27},totalExps{y,28}, 'm', 1);
%         p14a=shadedErrorBar(totalExps{z,26},totalExps{z,27},totalExps{z,28}, 'b', 1);
%         xlabel('time post-anaphase onset (min)') 
%         ylabel(strcat('\Delta distance between centroids (',micron,')'))
%         axis([xmin xlim ymin ylim])
%         legend([p13.patch,p14.patch,p13a.patch,p14a.patch],'shSPAST ctrl DMSO','shSPAST ctrl SPZ','shSPAST dox DMSO', 'shSPAST dox SPZ','Location','northwest')
%         legend('boxoff')
%         set(gca,'FontSize', 16);
% 
%         
%     else
%     end
%     
% savefig(testTitle);

% close all


