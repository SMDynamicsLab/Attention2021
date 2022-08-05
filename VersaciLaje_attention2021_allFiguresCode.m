clear; close all
%exp(year,attention).name(subject).   %%year=1: 2018; year=2: 2019; %attention=1: NORMAL; attention=2: HIGH
%                 .antena(subject)
%                 .workload(subject,feedback)  %%feedback=1: without feedback; feedback=2: with feedback
%                 .errors(subject).errorType{feedback}
%                 .datos(subject,feedback,pert).asy(trial,tap)  %%pert=1: -60; pert=2: +60
%                                             .resinchronizationArea
%                 .ERPs(subject,feedback,electrode).respLocked(epoch,sample)
%                                                 .originalStimLocked(epoch,sample)
%                                                 .respBaselines(trial)
%                                                 .stimBaselines(trial)


load VersaciLaje_attention2021_data
plotedElectrodes=[1:4 11:14]; %centro-frontal electrodes

% Set path to print figures to a file
path='C:\Users\focod\OneDrive\Escritorio\'; %pc personal windows
printOn=0; %print figures to file, 1=On; 0=Off


%% Generate rawRespERPs

yeari=2; %only 2019
for attentioni=1:2
    for subjecti=1:NSujetos
        for feedbacki=1:2            
            for electrodei=1:14
                for perti=1:2
                    for triali=1:NTrials
                        for epochi=fisrtCuttedResp:NEpochs_max
                            rawRespERPs_zeroing(attentioni,subjecti,feedbacki,perti,electrodei,triali,epochi,:)=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).respERPs(electrodei,triali,epochi,:);
                        end
                    end
                end
            end
        end
    end
end

%% Stim Invalid Epochs and Electrodes

q1=asyMin;
q2=asyMax;

%invalidEpochs_prePert(attentioni,subjecti,feedbacki,perti,electrode,triali,epochi)
%nan: without epoch. 0: Valid Epoch; 1: max 75 microV; 2: doesn't survive to Step Function algorithm; 

invalid_stim_temp=zeros(2,NSujetos,2,2,14,NTrials,NStim_prePert);
sample_freq=128; %sampling frequency

%Invalid epochs from begining and pre pert
invalidEpochs_prePert=invalidEpochs_fromBegining(:,:,:,:,:,:,end-NStim_prePert+1:end);
          
%Invalid electrodes pre pert
%invalidElectrodes_prePert(attention,NSujetos,feedback,electrode)
%0: valid electrode; 1: invalid in raw Resp Cut Process 
%2: bigger than maximun tolerable value; 3: less than minValidEpochsTolerable_perElectrode
invalidElectrodes_temp=zeros(2,NSujetos,2,14);
NEpoch=NTrials*NStim_prePert*2; %Epochs in pre pert zone

for attentioni=1:2
    for subjecti=1:NSujetos
        for feedbacki=1:2
            for electrodei=1:14
                counter=0;
                for perti=1:2
                    for triali=1:NTrials
                        for epochi=1:NStim_prePert
                            y=invalidEpochs_prePert(attentioni,subjecti,feedbacki,perti,electrodei,triali,epochi);

                            if y~=0
                                counter=counter+1;
                            end
                        end
                    end
                end
  
                if counter/NEpoch*100>maxTolerableTrialsRechazados
                    invalidElectrodes_temp(attentioni,subjecti,feedbacki,electrodei)=2;
                end
            end
        end
    end
end

%All invalid electrodes
invalidElectrodes_temp(find(fisrtCuttedResp==1))=1;
invalidElectrodes_prePert=invalidElectrodes_temp;

%% Cut Stim Locked Material
yeari=2;
for attentioni=1:2
    for subjecti=1:NSujetos
        for feedbacki=1:2
            for perti=1:2
                for electrodei=1:14
                    for triali=1:NTrials
                        i=0;
                        for epochi=NEpochs_max-NStim_prePert+1:NEpochs_max %recorta en zona pre pert
                            i=i+1;
                            asy=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asyPrePertFromBegining(triali,epochi);
                            %Verify that asy fall between percentiles
                            if asy>=asyMin && asy<=asyMax
                                ind1=-respCutWindow(1)+1+stimCutWindow(1)-asy;
                                ind2=-respCutWindow(1)+1+stimCutWindow(2)-asy;
                                recorte=squeeze(rawRespERPs_zeroing(attentioni,subjecti,feedbacki,perti,electrodei,triali,epochi-fisrtCuttedResp+1,ind1:ind2));
                                l=length(recorte);
                                %Decide if epoch is valid or not
                                y=invalidEpochs_prePert(attentioni,subjecti,feedbacki,perti,electrodei,triali,i);
                                if y~=0
                                    recorte=nan(1,l);
                                end
                            else
                                recorte=nan(1,l);
                            end
                            rawStimERPs(attentioni,subjecti,feedbacki,perti,electrodei,triali,i,:)=recorte;
                        end
                    end
                end
            end
        end
    end
end



%% Vector calculation

%Perturbations pooled
stimERPs_temp1(:,:,:,:,1:NTrials,:,:)=rawStimERPs(:,:,:,1,:,:,:,:);
stimERPs_temp1(:,:,:,:,NTrials+1:NTrials*2,:,:)=rawStimERPs(:,:,:,2,:,:,:,:);
L=size(stimERPs_temp1,7);

%Mean Trials and Epochs
stimERPs_temp2=squeeze(nanmean(reshape(stimERPs_temp1,[2,NSujetos,2,14,NTrials*2*NStim_prePert,L]),5));

%Discard invalid electrodes
for attentioni=1:2
    for subjecti=1:NSujetos
        for feedbacki=1:2
            for electrodei=1:14
                y=invalidElectrodes_prePert(attentioni,subjecti,feedbacki,electrodei);
                if y~=0 %invalid electrode
                    stimERPs_temp2(attentioni,subjecti,feedbacki,electrodei,:)=nan;
                end
            end
        end
    end
end

%Average electrtode         
stimERPs_temp3=squeeze(nanmean(stimERPs_temp2(:,:,:,plotedElectrodes,:),4));

% AdJar
%Perturbations pooled
respERPs_temp1(:,:,:,:,1:NTrials,:,:)=rawRespERPs_zeroing(:,:,:,1,:,:,:,:);
respERPs_temp1(:,:,:,:,NTrials+1:NTrials*2,:,:)=rawRespERPs_zeroing(:,:,:,2,:,:,:,:);

%Mean Trials and Epochs
NSamples=respCutWindow(2)-respCutWindow(1)+1;
respERPs_temp2=squeeze(nanmean(reshape(respERPs_temp1,[2,NSujetos,2,14,NTrials*2*NEpochs_max,NSamples]),5));

%Name unification
rERPs_adjar=respERPs_temp2;
sERPs_adjar=stimERPs_temp2;

%NAsy calculation
yeari=2;
for attentioni=1:2
    for subjecti=1:NSujetos
        for feedbacki=1:2
            for perti=1:2
                for electrodei=1:14
                    for triali=1:NTrials
                        i=0;
                        for epochi=NEpochs_max-NStim_prePert+1:NEpochs_max %recorta en zona pre pert
                            i=i+1;
                            temp=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asyPrePertFromBegining(triali,epochi);
                            asy_temp(attentioni,subjecti,feedbacki,perti,electrodei,triali,i)=temp;
                        end
                    end
                end
            end
        end
    end
end
asy_temp(find(invalidEpochs_prePert~=0))=nan;

asy=[];
diffStim=[];
for attentioni=1:2
    for subjecti=1:NSujetos
        for feedbacki=1:2
            for electrodei=1:14
                %Calculate NAsy for each electrode (max=NStim_prePert*NTrials_correctosPerBlock)
                asy=squeeze(asy_temp(attentioni,subjecti,feedbacki,:,electrodei,:,:));
                asy=asy(:);
                asy=asy(~isnan(asy)); 
                NAsy_actual=length(asy);%Substract Resp Locked Material
                NAsy(attentioni,subjecti,feedbacki,electrodei)=NAsy_actual;
                newStim=squeeze(sERPs_adjar(attentioni,subjecti,feedbacki,electrodei,:));
                diffStim=zeros(1,length(newStim));
                resp=squeeze(rERPs_adjar(attentioni,subjecti,feedbacki,electrodei,:))/NAsy_actual;
                for i=1:length(asy) %scan all asynchronies
                    ind1=abs(respCutWindow(1))+1-asy(i)+stimCutWindow(1);
                    ind2=abs(respCutWindow(1))+1-asy(i)+stimCutWindow(2);
                    recorteRespAlineado=resp(ind1:ind2);
                    newStim=[newStim-recorteRespAlineado];
                    diffStim=[diffStim+recorteRespAlineado'];
                end

                    stimERPs_corrected(electrodei,:,feedbacki,subjecti,attentioni)=newStim;
                    diffBetweenStimAndStim_corrected(electrodei,:,feedbacki,subjecti,attentioni)=diffStim;
           end
        end
    end
end

%Diff ERP 
NSamples_stimCutWindow=stimCutWindow(2)-stimCutWindow(1)+1;
diffERP=permute(squeeze(nanmean(diffBetweenStimAndStim_corrected(plotedElectrodes,1:NSamples_stimCutWindow,:,:,:))),[3 1 2 4]);

%Corrected Stim Locked ERP 
correctedERPs=permute(squeeze(nanmean(stimERPs_corrected(plotedElectrodes,1:NSamples_stimCutWindow,:,:,:))),[3 1 2 4]);
leftBorder=stimCutWindow(1);
rightBorder=stimCutWindow(2);

% Original ERPs 
%Average Feedback Conditions
originalERPs2=permute(squeeze(mean(stimERPs_temp3,3)),[2 3 1]);

%Corrected ERPs
%Average of feedback conditions
correctedERPs2=squeeze(mean(correctedERPs,3));

%Diff ERPS
%Average feedback conditions
diffERP_fbkAveraged=squeeze(mean(diffERP,3));

%Resp ERP
%Pool perturbations
temp=[];
temp(:,:,:,:,1:NTrials,:,:)=rawRespERPs_zeroing(:,:,:,1,:,:,:,:);
temp(:,:,:,:,NTrials+1:NTrials*2,:,:)=rawRespERPs_zeroing(:,:,:,2,:,:,:,:);
%Mean epochs, trials
temp2=squeeze(nanmean(reshape(temp,[2 NSujetos 2 14 NTrials*2*NEpochs_max NSamples]),5));
%Mean ploted electrodes
temp3=squeeze(nanmean(temp2(:,:,:,plotedElectrodes,:),4));
% respERP=squeeze(nanmean(nanmean(nanmean(temp(:,:,:,plotedElectrodes,:,:,:),6),5),4));
respERP=temp3;
%Average feedback conditions
respERP_fbkAveraged=squeeze(mean(respERP,3));
%Average attention conditions
% respERP_fbkAveraged=permute(squeeze(mean(respERP,1)),[2 1 3]);


%% Temporal series and MA arrays 
pre=[];
pos=[];
T=[];
M=[];

for yeari=1:2
    for attentioni=1:2
        for subjecti=1:NSujetos
            for feedbacki=1:2
                for perti=1:2
                    asy_actual=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy(:,prePertZone);
                    pre(yeari,attentioni,subjecti,feedbacki,perti)=mean(mean(asy_actual,2));
                    pos(yeari,attentioni,subjecti,feedbacki,perti)=mean(mean(exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy(:,posPertZone),2)); 
                    preSTD(yeari,attentioni,subjecti,feedbacki,perti)=mean(std(exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy(:,prePertZone),0,2)); 
                    posSTD(yeari,attentioni,subjecti,feedbacki,perti)=mean(std(exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy(:,posPertZone),0,2)); 
                    series(yeari,attentioni,subjecti,feedbacki,perti,:)=nanmean(exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy);
                    series_ES(yeari,attentioni,subjecti,feedbacki,perti,:)=nanstd(exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy)/sqrt(NTrials);
                    STDInter_temp(yeari,attentioni,subjecti,feedbacki,perti)=std(nanmean(exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy(:,prePertZone),2));
                    STDIntra_temp(yeari,attentioni,subjecti,feedbacki,perti)=mean(nanstd(exp(yeari,attentioni).datos(subjecti,feedbacki,perti).asy(:,prePertZone),0,2));
                    accumulatedAsynchrony_temp(yeari,attentioni,subjecti,feedbacki,perti)=exp(yeari,attentioni).datos(subjecti,feedbacki,perti).resynchronizationArea;
                end
            end
        end
    end
end


%% Color Setup 

NColors=48;
indColors1=[1 13 27 35]; %posibilidad5
brilloInicial=0.4;
col_temp=ametrine(NColors);
col=col_temp(indColors1,:);

%% FIGURE 5 - Auditory ERPs and workload

h5=figure(5); 
clf(5);
lwidth=1;
msize=3;
fsize = 10;
fig_size_cm=[20 7];

set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

subplot(10,10,[1:6 11:16 21:26 31:36 41:46 51:56 61:66 71:76 81:86]) 

ptp=[];

%Specific Colors setup
deltaBrillo=0.5;
indColors2=[floor(mean([indColors1(1) indColors1(2)])); floor(mean([indColors1(3) indColors1(4)]))];

brillo=brilloInicial;
for i=1:2
    col1=brighten(ametrine(NColors),-brillo);
    col2(i,:)=col1(indColors2(i),:);
    brillo=brillo-deltaBrillo; 
end


%Significance Zone - 1 tale
x0=133;
xf=156;

transp=0.25;
col=[1 1 1]*0.4;
patch([x0 x0 xf xf],[-5 5 5 -5],col,'facealpha',transp,'edgecolor','none','edgealpha',transp); hold on
ylabel(['Voltage (' char(181) 'V)' ],'fontsize',fsize)
xlabel(['Time (ms)' ],'fontsize',fsize,'position',[130 -2.7 0])


%Plot
x=[leftBorder:1:rightBorder]*1000/sample_freq;
plot([stimCutWindow(1)*1000/128 stimCutWindow(2)*1000/128],[0 0],'k-'); hold on
plot([0 0],[-20 20],'k-'); hold on

for attentioni=1:2
    y=mean(correctedERPs2(:,:,attentioni));
    ES=std(correctedERPs2(:,:,attentioni))/sqrt(NSujetos);
    ptp(attentioni)=plot([leftBorder:1:rightBorder]*1000/128,y,'-','color',[col2(attentioni,:)],'linewidth',lwidth); hold on
    lineProps.width = 1;
    lineProps.edgestyle = '-';
    lineProps.col{1}=[col2(attentioni,:)];
    mseb(x,y,ES,lineProps,0.8);
end
legend(ptp,{'NORMAL','HIGH'},'position',[0.48 0.28 0.01 0.01],'fontsize',fsize-1)
legend boxoff
set(gca,'ytick',[-3 -2 -1 0 1 2 3],'fontsize',fsize)

ylim([-2.2 1.8])
xlim([stimCutWindow(1)*1000/128 stimCutWindow(2)*1000/128])
set(gca,'xtick',[-100 0 100 200 300],'fontsize',fsize)
set(gca,'box','off')
set(gca,'xcolor','k')
text(-0.11,0.975,'A','units','normalized','fontsize',fsize+3)

% Workload
subplot(10,10,[8:10 18:20 28:30 38:40 48:50 58:60 68:70 78:80 88:90 98:100 ]) 

yeari=2;
workload=nan(NSujetos,2,2);
for attentioni=1:2
    workload(:,:,attentioni)=exp(yeari,attentioni).workload;
end

%Mean
y=[mean(workload(:,:,1),2) mean(workload(:,:,2),2)];

%Boxplot - Feedback Conditions Averaged
ptb = boxplot(y,'colors',col2);
set(gca,'xtick',[1:2],'xticklabel',{'NORMAL','HIGH'},'fontsize',fsize)
set(gca,'ytick',[0 2 4 6 8],'fontsize',fsize)
ylim([1 9])
ylabel('Likert Scale (1-9)','fontsize',fsize) %,'position',[0.3 4.6 0])
p = kruskalwallis(y,[],'off');    
set(findobj(gca,'type','line'),'linew',1)
set(gca,'box','off')
text(-0.21,0.975,'B','units','normalized','fontsize',fsize+3)
set(ptb,'linestyle','-')
set(ptb(7,1),'MarkerEdgeColor',col2(1,:))

%Scatter plot
%Manually set horizontal position of points to avoid superposition
d=0.035;
ypos(:,1)=[0 d 0 0 d 0 -d -d 0 d -d];
ypos(:,2)=[d 0 0 -d d -d 0 d d -d -d];
for barrai=1:2
    hold on
    for i=1:NSujetos
        hhh=plot(barrai+0.35+ypos(i,barrai),y(i,barrai),'o','markeredgecolor','k','markersize',2.5,'linewidth',1);
        drawnow;
        hhh.MarkerHandle.LineWidth = 0.1;
    end
end

%%Sifnificance bracket - Workload
yy=[8]; %bracket y coordinate 
s=12; %font size
w=0.5; %line width
for i=1:length(yy)
    orientacion=0;  %0 down,1 up
    l=1; %bracket length
    x0=(i-1)*l+1; %initial position in x
    y0=yy(i); %initial position in y
    yl=ylim; %y limits
    b=(yl(2)-yl(1))/45; %small lines length
    linea=[repmat('_',1,l*31)];
    if orientacion==1
        b=-b;
    end
    if mod(i,2)==1
        texto='*';
    else
        texto='HIGH';
    end    
    line([x0 x0+l],[y0 y0],'linewidth',w,'color','k')
    line([x0 x0],[y0 y0-b],'linewidth',w,'color','k')
    line([x0+l x0+l],[y0 y0-b],'linewidth',w,'color','k')
    text(x0+l/2+0,y0+b,texto,'fontsize',s)    
end


if printOn==1
    print(h5,'-depsc','-painters',[path 'figERPsWorkload.eps']);
end

return

%% FIGURE A2 (appendix) - Resp locked ERPs, Original ERPS and Corrected ERPs 

h10=figure(10); 
clf(10);
fig_size_cm=[15 6];
fsize = 9;
lwidth=1;
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

% Plot Stim ERPs
subplot(10,3,[1 4 7 10 13 16 19 22 25]+1)
ptp=[];

%Average of feedback conditions
x=[leftBorder:1:rightBorder]*1000/128;
line([stimCutWindow(1)*1000/128 stimCutWindow(2)*1000/128],[0 0],'linestyle','-','color','k')
line([0 0],[-20 20],'linestyle','-','color','k')
hold on
for attentioni=1:2
    y=mean(originalERPs2(:,:,attentioni));
    ES=std(originalERPs2(:,:,attentioni))/sqrt(NSujetos);
    ptp(attentioni)=plot([leftBorder:1:rightBorder]*1000/128,y,'-','color',[col2(attentioni,:)],'linewidth',lwidth); hold on
    lineProps.width = 1;
    lineProps.edgestyle = '-';
    lineProps.col{1}=[col2(attentioni,:)];
    mseb(x,y,ES,lineProps,0.8);
end
ylim([-3.5 3.4])
xlim([stimCutWindow(1)*1000/128 stimCutWindow(2)*1000/128])
set(gca,'ytick',[])
set(gca,'xtick',[])
set(gca,'xtick',[-100 0 100 200 300],'fontsize',fsize)
set(gca,'ytick',[-3 -2 -1 0 1 2 3],'yticklabel',[])
set(gca,'box','off')
title(['Stim-locked ERPs'],'fontsize',fsize-1,'fontweight','normal')
text(-0.18,1.06,'B','fontsize',fsize+1,'units','normalized')


%Corrected ERPs
subplot(10,3,[1 4 7 10 13 16 19 22 25]+2)
ptp=[];

%Significance Zone - 1 tale
x0=133;
xf=156;
transp=0.25;
col=[1 1 1]*0.4;
patch([x0 x0 xf xf],[-5 5 5 -5],col,'facealpha',transp,'edgecolor','none','edgealpha',transp)
hold on
line([stimCutWindow(1)*1000/128 stimCutWindow(2)*1000/128],[0 0],'linestyle','-','color','k')
line([0 0],[-20 20],'linestyle','-','color','k')
hold on
for attentioni=1:2
    y=mean(correctedERPs2(:,:,attentioni));
    ES=std(correctedERPs2(:,:,attentioni))/sqrt(NSujetos);
    ptp(attentioni)=plot([leftBorder:1:rightBorder]*1000/128,y,'-','color',[col2(attentioni,:)],'linewidth',lwidth); hold on
    lineProps.width = 1;
    lineProps.edgestyle = '-';
    lineProps.col{1}=[col2(attentioni,:)];
    mseb(x,y,ES,lineProps,0.8);
end
set(gca,'ytick',[-3 -2 -1 0 1 2 3],'yticklabel',[])
ylim([-3.5 3.4])
xlim([stimCutWindow(1)*1000/128 stimCutWindow(2)*1000/128])
set(gca,'xtick',[-100 0 100 200 300],'fontsize',fsize)
set(gca,'box','off')
text(-0.18,1.06,'C','fontsize',fsize+1,'units','normalized')
title(['Corrected Stim-locked ERPs'],'fontsize',fsize-1,'fontweight','normal','position',[200 3.47 0])
leg=legend(ptp,{'NORMAL','HIGH'},'position',[0.89 0.270 0.01 0.01]);
leg.ItemTokenSize = [15,20];
legend boxoff

%Resp ERPs - Feedback conditions averaged
subplot(10,3,[1 4 7 10 13 16 19 22 25]+0)
leftBorder_resp=respCutWindow(1);
rightBorder_resp=respCutWindow(2);
x=[leftBorder_resp:1:rightBorder_resp]*1000/sample_freq;
line([respCutWindow(1)*1000/sample_freq respCutWindow(2)*1000/sample_freq],[0 0],'linestyle','-','color','k')
line([0 0],[-10 10],'linestyle','-','color','k')
hold on
for attentioni=1:2
    y=squeeze(mean(respERP_fbkAveraged(attentioni,:,:),2));
    ES=squeeze(std(respERP_fbkAveraged(attentioni,:,:),0,2))/sqrt(NSujetos);
    ptp(attentioni)=plot([leftBorder_resp:1:rightBorder_resp]*1000/128,y,'-','color',[col2(attentioni,:)],'linewidth',lwidth); hold on
    lineProps.width = 1;
    lineProps.edgestyle = '-';
    lineProps.col{1}=[col2(attentioni,:)];
    mseb(x,y',ES',lineProps,0.8);
end
xlim([-150 350])
ylim([-3.5 3.4])
set(gca,'xtick',[-100 0 100 200 300],'fontsize',fsize)
set(gca,'ytick',[-3 -2 -1 0 1 2 3],'fontsize',fsize)
set(gca,'box','off')
title(['Resp-locked ERPs'],'fontsize',fsize-1,'fontweight','normal')
text(-0.32,1.06,'A','fontsize',fsize+1,'units','normalized')

%General Y axis
set(gcf,'NextPlot','add');
axes;
h = ylabel(['Voltage (' char(181) 'V)' ],'fontsize',fsize);
set(gca,'Visible','off');
set(h,'Visible','on');

%%General X axis
set(gcf,'NextPlot','add');
axes;
h = xlabel(['Time (ms)'],'fontsize',fsize,'position',[0.5 0 0]);
set(gca,'Visible','off');
set(h,'Visible','on');

set(gca,'LooseInset',get(gca,'TightInset'));

if printOn==1
    print(h10,'-depsc','-painters',[path 'figresultadoTapAdjar.eps']);
end


%% FIGURE 2 - Temporal Series of asynchronies - MA and SD 

deltaBrillo=0.3; %variación de brillo entre una condición y la siguiente
h1=figure(1);
clf(1);
lwidth=1;
msize=10;
fsize = 9;
fsizelegend = 8;
fig_size_cm=[20 8];
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

% SUBPLOT - Temporal series
o=[1:4 13:16 25:28];
indPlot=[o; o+36; o+72; o+108; o+144];

ploti=0;
for attentioni=1:2
    for feedbacki=1:2
        ploti=ploti+1;
        subplot(12,12,indPlot(ploti,:))

        line([0 NTaps+1],[0 0],'LineStyle','-','Color','k'); hold on
            line([8 8],[-220 200],'LineStyle','-','Color','k');

        for perti=1:2
            if perti==1
                linea='-';
                if attentioni==2 && feedbacki==2
                brillo=brilloInicial-0.25;
                else
                    brillo=brilloInicial-0.45;
                end
                for i=1:4
                    col1=brighten(ametrine(NColors),-brillo);
                    col(i,:)=col1(indColors1(i),:);
                    brillo=brillo-deltaBrillo;
                end
            else
                linea='-';
                brillo=brilloInicial;
                for i=1:4
                    col1=brighten(ametrine(NColors),-brillo);
                    col(i,:)=col1(indColors1(i),:);
                    brillo=brillo-deltaBrillo;
                end
            end
            for yeari=1:2
                for subjecti=1:NSujetos
                    asy=squeeze(series(yeari,attentioni,subjecti,feedbacki,perti,:));
                    ES=squeeze(series_ES(yeari,attentioni,subjecti,feedbacki,perti,:));
                    ptp(ploti)=plot(1:NTaps,asy,'.','linestyle',linea,'MarkerSize',msize,'LineWidth',lwidth,'MarkerFaceColor',[col(ploti,:)],'Color',[col(ploti,:)]);
                    hold on
                end
            end
            
            %Axis limits
            xlim([0 NTaps+1])
            ylim([-220 110])
            %Ticks
            set(gca, 'YTick',[-100 0 100]);
            set(gca,'box','off') 
            set(gca,'xtick',[])
            set(gca,'XTicklabel',[])
            set(gca,'xcolor', 'w')
            set(gca,'ycolor', 'k')
            set(gca,'FontSize', fsize);
            
            %Conditions Labels
            if ploti==1
                texto='NORMAL - No feedback';
                text(-0.28,0.99,'A','units','normalized','fontsize',fsize)
            elseif ploti==2
                texto='NORMAL - With feedback';
            elseif ploti==3
                texto='HIGH - No feedback';
            else
                texto='HIGH - With feedback';
            end
                      
        end
   end
end
set(gca,'xtick',[3 8 13 18],'xticklabel',[-5 0 5 10],'xcolor','k','fontsize',fsize)

% Legend
brillo=brilloInicial;
for i=1:4
    col1=brighten(ametrine(NColors),-brillo);
    col(i,:)=col1(indColors1(i),:);
    brillo=brillo-deltaBrillo;
end
for ploti=1:4
    ptp(ploti)=plot(100,0,'o','linestyle',linea,'MarkerSize',msize/2,'LineWidth',2*lwidth,'MarkerFaceColor',[col(ploti,:)],'Color',[col(ploti,:)]);
end
legend(ptp,{'NORMAL, No FBK','NORMAL, With FBK','HIGH, No FBK','HIGH, With FBK'},'position',[0.62 0.1 0.1 0.01],'fontsize',fsizelegend) %[0.7 0.85 0.1 0.01] referencia arriba
legend boxoff

% Subplot MA - MA (year, att, subj, fbk, pert)
T1=squeeze(pre(1,1,:,1,1));
T2=squeeze(pre(1,1,:,1,2));
T3=squeeze(pre(1,1,:,2,1));
T4=squeeze(pre(1,1,:,2,2));
T5=squeeze(pre(1,2,:,1,1));
T6=squeeze(pre(1,2,:,1,2));
T7=squeeze(pre(1,2,:,2,1));
T8=squeeze(pre(1,2,:,2,2));
T9=squeeze(pre(2,1,:,1,1));
T10=squeeze(pre(2,1,:,1,2));
T11=squeeze(pre(2,1,:,2,1));
T12=squeeze(pre(2,1,:,2,2));
T13=squeeze(pre(2,2,:,1,1));
T14=squeeze(pre(2,2,:,1,2));
T15=squeeze(pre(2,2,:,2,1));
T16=squeeze(pre(2,2,:,2,2));
%Mean Perturbation Conditions
TT1=mean([T1,T2],2);
TT2=mean([T3,T4],2);
TT3=mean([T5,T6],2);
TT4=mean([T7,T8],2);
TT5=mean([T9,T10],2);
TT6=mean([T11,T12],2);
TT7=mean([T13,T14],2);
TT8=mean([T15,T16],2);
%Pool Year
NMA_pre=[[TT1;TT5],[TT2;TT6],[TT3;TT7],[TT4;TT8]];

%Specific color setup
brillo=brilloInicial;
for i=1:4
    col1=brighten(ametrine(NColors),-brillo);
    col(i,:)=col1(indColors1(i),:);
    brillo=brillo-deltaBrillo; 
end

%Bar Plot
o=[6:8 18:20 30:32 42:44 54:56 66:68 78:80 90:92 102:104]+12; %Paneles B y C en columnas y más estirados
subplot(12,12,o)
y=mean(NMA_pre);
d=0.06; %distance betwenn groups
xbar=[1+d 2-d 3+d 4-d ];
ES=std(NMA_pre)/sqrt(2*NSujetos);
for i=1:4
    bar(xbar(i),y(i),'facecolor',col(i,:)); hold on
    errorbar(xbar(i),y(i),ES(i),'-','color',col(i,:),'linewidth',lwidth)
end

%Scatter plot
for barrai=1:4
    hold on
    for i=1:NSujetos*2
        plot(xbar(barrai)+0.25,NMA_pre(i,barrai),'o','markerfacecolor','w','markersize',2.5,'linewidth',1,'markeredgecolor','k');
        hhh=plot(xbar(barrai)+0.25,NMA_pre(i,barrai),'o','markeredgecolor','k','markersize',2.5,'linewidth',1);
        drawnow;
        hhh.MarkerHandle.LineWidth = 0.05;
    end
end

xlim([0 5])
ylim([-135 45]) 
set(gca,'ytick',[-100 -50 0 50],'xtick',[],'fontsize',fsize)
title('\it MA_P_R_E','fontsize',fsize,'Interpreter','tex','fontweight','normal')
ylabel('(ms)','fontsize',fsize)
text(-0.3,1.09,'B','units','normalized','fontsize',fsize)
                
%Significance bracket
y=[30]; %bracket y coordinate
s=15; %font size
w=1; %line width
for i=1:length(y)
    orientacion=0;  %0 down,1 up
    l=2; %length
    x0=(i-1)*l+1.5; %initial position x
    y0=y(i); %initial position x
    yl=ylim; 
    b=(yl(2)-yl(1))/45; %small lines length
    linea=[repmat('_',1,l*31)];
    if orientacion==1
        b=-b;
    end
    if mod(i,2)==1
        texto='**';
    else
        texto='HIGH';
    end    
    line([x0 x0+l],[y0 y0],'linewidth',w,'color','k')
    line([x0 x0],[y0 y0-b],'linewidth',w,'color','k')
    line([x0+l x0+l],[y0 y0-b],'linewidth',w,'color','k')
    text(x0+l/2-0.2,y0+2,texto,'fontsize',s)    
end

%ANOVA
factorA=[ones(NSujetos*2,1);2*ones(NSujetos*2,1)];
Bnames={'FEEDBACK'}; %factor B
Anames={'ATENCION'}; %factor A
M=[NMA_pre(:,1:2);NMA_pre(:,3:4)];
[tbl,rm]=simple_mixed_anova(M,factorA,Bnames,Anames);


% SUBPLOT - SD 
o=[6:8 18:20 30:32 42:44 54:56 66:68 78:80 90:92 102:104]+4+12; 
subplot(12,12,o)

temp=STDIntra_temp;
T1=squeeze(temp(1,1,:,1,1));
T2=squeeze(temp(1,1,:,1,2));
T3=squeeze(temp(1,1,:,2,1));
T4=squeeze(temp(1,1,:,2,2));
T5=squeeze(temp(1,2,:,1,1));
T6=squeeze(temp(1,2,:,1,2));
T7=squeeze(temp(1,2,:,2,1));
T8=squeeze(temp(1,2,:,2,2));
T9=squeeze(temp(2,1,:,1,1));
T10=squeeze(temp(2,1,:,1,2));
T11=squeeze(temp(2,1,:,2,1));
T12=squeeze(temp(2,1,:,2,2));
T13=squeeze(temp(2,2,:,1,1));
T14=squeeze(temp(2,2,:,1,2));
T15=squeeze(temp(2,2,:,2,1));
T16=squeeze(temp(2,2,:,2,2));
%Mean Perturbation Conditions
TT1=mean([T1,T2],2);
TT2=mean([T3,T4],2);
TT3=mean([T5,T6],2);
TT4=mean([T7,T8],2);
TT5=mean([T9,T10],2);
TT6=mean([T11,T12],2);
TT7=mean([T13,T14],2);
TT8=mean([T15,T16],2);
%Pool Year
STDIntra=[[TT1;TT5],[TT2;TT6],[TT3;TT7],[TT4;TT8]];
B=STDIntra;

y=mean(B);
ES=std(B)/sqrt(2*NSujetos);
for i=1:4
    bar(xbar(i),y(i),'facecolor',col(i,:)); hold on
    errorbar(xbar(i),y(i),ES(i),'-','color',col(i,:),'linewidth',lwidth)
end

%Scatter plot
for barrai=1:4
    hold on
    for i=1:NSujetos*2
        plot(xbar(barrai)+0.25,STDIntra(i,barrai),'o','markerfacecolor','w','markersize',2.5,'linewidth',1,'markeredgecolor','k');
        hhh=plot(xbar(barrai)+0.25,STDIntra(i,barrai),'o','markeredgecolor','k','markersize',2.5,'linewidth',1);
        drawnow;
        hhh.MarkerHandle.LineWidth = 0.1;
    end
end
xlim([0 5])
ylim([0 30]) 
set(gca,'ytick',[0 10 20 30],'xtick',[],'fontsize',fsize)
title('\it SD_P_R_E','fontsize',fsize,'Interpreter','tex','position',[2.45 30.5 0],'fontweight','normal')

ylabel('(ms)','fontsize',fsize)
text(-0.23,1.09,'C','units','normalized','fontsize',fsize)

%ANOVA
% factorA=[ones(NSujetos*2,1);2*ones(NSujetos*2,1)];
% Bnames={'FEEDBACK'}; %factor B
% Anames={'ATENCION'}; %factor A
% M=[B(:,1:2);B(:,3:4)];
% [tbl,rm]=simple_mixed_anova(M,factorA,Bnames,Anames);

% text(0.65,0.35,['feed p=' num2str(round(1000*tbl.pValue(4))/1000) ],'Units','normalized','Color','black','FontSize',10)
% text(0.65,0.25,['ate p=' num2str(round(100000*tbl.pValue(2))/100000) ],'Units','normalized','Color','black','FontSize',10)
% text(0.65,0.15,['Feed X Ate p=' num2str(round(1000*tbl.pValue(5))/1000) ],'Units','normalized','Color','black','FontSize',10)

% General Titles and Axes
%Yaxe
set(gcf,'NextPlot','add');
axes;
h = ylabel(['Asynchrony {\it e_n} (ms)'] ,'fontsize',fsize,'position',[-0.06 0.5 0],'Interpreter','tex');
set(gca,'Visible','off');
set(h,'Visible','on');

%Xaxe
set(gcf,'NextPlot','add');
axes;
h = xlabel(['Beep number{\it n} (relative to perturbation beep)'] ,'fontsize',fsize,'position',[0.13 -0.065 0],'Interpreter','tex');
set(gca,'Visible','off');
set(h,'Visible','on');

if printOn==1
    print(h1,'-depsc','-painters',[path 'figure1.eps']);
end

return


%% FIGURE 3 - Resynchronization efficiency
h2=figure(2);
clf(2);
subplot(11,10,[1:9 11:19 21:29 31:39 41:49 51:59 61:69 71:79 81:89 91:99])
lwidth=1;
msize=3;
fsize = 8;
fig_size_cm=[10 5];
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

%Specific color setup
indColors1=[1 1 13 13 27 27 35 35]; 
brillo=brilloInicial;
for i=[1 3 5 7]
    col1=brighten(ametrine(NColors),-brillo);
    col(i,:)=col1(indColors1(i),:);
    col(i+1,:)=col1(indColors1(i),:);
    brillo=brillo-deltaBrillo; 
end

%Accumulated asynchrony 
T1=squeeze(accumulatedAsynchrony_temp(1,1,:,1,1));
T2=squeeze(accumulatedAsynchrony_temp(1,1,:,1,2));
T3=squeeze(accumulatedAsynchrony_temp(1,1,:,2,1));
T4=squeeze(accumulatedAsynchrony_temp(1,1,:,2,2));
T5=squeeze(accumulatedAsynchrony_temp(1,2,:,1,1));
T6=squeeze(accumulatedAsynchrony_temp(1,2,:,1,2));
T7=squeeze(accumulatedAsynchrony_temp(1,2,:,2,1));
T8=squeeze(accumulatedAsynchrony_temp(1,2,:,2,2));
T9=squeeze(accumulatedAsynchrony_temp(2,1,:,1,1));
T10=squeeze(accumulatedAsynchrony_temp(2,1,:,1,2));
T11=squeeze(accumulatedAsynchrony_temp(2,1,:,2,1));
T12=squeeze(accumulatedAsynchrony_temp(2,1,:,2,2));
T13=squeeze(accumulatedAsynchrony_temp(2,2,:,1,1));
T14=squeeze(accumulatedAsynchrony_temp(2,2,:,1,2));
T15=squeeze(accumulatedAsynchrony_temp(2,2,:,2,1));
T16=squeeze(accumulatedAsynchrony_temp(2,2,:,2,2));
accumulatedAsynchrony=[[T1;T9],[T2;T10],[T3;T11],[T4;T12],[T5;T13],[T6;T14],[T7;T15],[T8;T16]];

%(asy perturb bip - post baseline)
asy8=series(:,:,:,:,:,8); %asy perturb bip
despla=asy8-pos;
T1=squeeze(despla(1,1,:,1,1));
T2=squeeze(despla(1,1,:,1,2));
T3=squeeze(despla(1,1,:,2,1));
T4=squeeze(despla(1,1,:,2,2));
T5=squeeze(despla(1,2,:,1,1));
T6=squeeze(despla(1,2,:,1,2));
T7=squeeze(despla(1,2,:,2,1));
T8=squeeze(despla(1,2,:,2,2));
T9=squeeze(despla(2,1,:,1,1));
T10=squeeze(despla(2,1,:,1,2));
T11=squeeze(despla(2,1,:,2,1));
T12=squeeze(despla(2,1,:,2,2));
T13=squeeze(despla(2,2,:,1,1));
T14=squeeze(despla(2,2,:,1,2));
T15=squeeze(despla(2,2,:,2,1));
T16=squeeze(despla(2,2,:,2,2));
desplazamiento=[[T1;T9],[T2;T10],[T3;T11],[T4;T12],[T5;T13],[T6;T14],[T7;T15],[T8;T16]];

%efficiency (e=(asy perturb bip - post baseline)/accumulatedAsynchrony)
e=desplazamiento./accumulatedAsynchrony;

% Bars
T=abs(e);
d=0.2;
xbar=[1-d 2-d 3-d 4-d 5+d 6+d 7+d 8+d];
q=sqrt(NSujetos*2);
for barrai=1:8
    ptp(barrai)=bar(xbar(barrai),mean(T(:,barrai)),'facecolor',col(barrai,:));
    hold on
    errorbar(xbar(barrai),mean(T(:,barrai)),std(T(:,barrai))/q,'-','Color',col(barrai,:),'linewidth',lwidth)
end

%Scatter plot
for barrai=1:8
    hold on
    for i=1:NSujetos*2
        plot(xbar(barrai)+0.25,abs(e(i,barrai)),'o','markerfacecolor','w','markersize',2.5,'linewidth',1,'markeredgecolor','k');
        hhh=plot(xbar(barrai)+0.25,abs(e(i,barrai)),'o','markeredgecolor','k','markersize',2.5,'linewidth',1);
        drawnow;
        hhh.MarkerHandle.LineWidth = 0.1;
    end
end
xlim([0 9]) 
ylim([0 1.05]) 
set(gca,'xtick',[],'xticklabel',[])
set(gca,'box','off','ytick',[0:0.1:1],'yticklabel',{'0','','','','','0.5','','','','','1'},'fontsize',fsize)
ylabel({'Resynchronization'; 'efficiency (adim)'},'fontsize',fsize)
for i=1:8
    if mod(i,2)==0
        texto='Pos';
    else
        texto='Neg';
    end
    text(xbar(i)-0.25,-0.06,texto,'fontsize',fsize-1)
end

%vector for ANOVA 4
M=[T(:,1);T(:,2);T(:,3);T(:,4);T(:,5);T(:,6);T(:,7);T(:,8)];
sujS=[repmat((1:NSujetos*2)',4,1);repmat((1:NSujetos*2)'+NSujetos*2,4,1)];
atencionA=[repmat(1,1,NSujetos*8)';repmat(2,1,NSujetos*8)'];
feedbackB=[repmat(1,1,NSujetos*4)';repmat(2,1,NSujetos*4)';repmat(1,1,NSujetos*4)';repmat(2,1,NSujetos*4)'];
pertD=repmat([repmat(1,1,NSujetos*2)';repmat(2,1,NSujetos*2)'],4,1);
nesting=[0 0 0 0 ; 0 0 0 0 ; 1 0 0 0 ; 0 0 0 0 ]; %SujS nested with atenciÃ³nA - pÃ¡gina 512, Keppel
[p2,table2,stats2,terms2]=anovan(M,{atencionA feedbackB sujS pertD},'random',3,'nested',nesting,'varnames',{'Aten' 'Feed' 'Sujs' 'Pert'},'model','full','display','off');

%Texts p values
% text(0.8,0.98,['Aten p=' num2str(round(10000*table2{2,7})/10000)],'units','normalized','fontsize',13)
% text(0.8,0.95,['Feed p=' num2str(round(10000*table2{3,7})/10000)],'units','normalized','fontsize',13)
% text(0.8,0.92,['Pert p=' num2str(round(10000*table2{5,7})/10000)],'units','normalized','fontsize',13)
% text(0.8,0.89,['Aten*Feed p=' num2str(round(10000*table2{6,7})/10000)],'units','normalized','fontsize',13)
% text(0.8,0.86,['Aten*Pert p=' num2str(round(10000*table2{7,7})/10000)],'units','normalized','fontsize',13)
% text(0.8,0.83,['Feed*Pert p=' num2str(round(10000*table2{9,7})/10000)],'units','normalized','fontsize',13)
% text(0.8,0.8,['Aten*Feed*Pert p=' num2str(round(10000*table2{11,7})/10000)],'units','normalized','fontsize',13)
% text(0.8,0.7,['NOTA: no se exluyÃ³ trial outlier'],'units','normalized','fontsize',13)
% text(0.8,0.67,['Ver resynchronizationEfficiency.m'],'units','normalized','fontsize',13)

%Significance factors and asteriks - on the right
text(1.0,0.55,'attention','fontsize',fsize-1,'units','normalized')
text(1.15,0.53,'*','fontsize',10,'units','normalized')
text(1.0,0.48,'feedback','fontsize',fsize-1,'units','normalized')
text(1.15,0.46,'**','fontsize',10,'units','normalized')
text(1.0,0.41,'pert. sign','fontsize',fsize-1,'units','normalized')
text(1.15,0.39,'*','fontsize',10,'units','normalized')

%Feedback texts
x01=(xbar(1)+xbar(2))/2-0.5;
x02=(xbar(3)+xbar(4))/2-0.6;
x03=(xbar(5)+xbar(6))/2-0.5;
x04=(xbar(7)+xbar(8))/2-0.6;
x=[x01 x02 x03 x04];
y0=-0.13;
fs=6;
text(x01, y0,'No FBK','fontsize',fs)
text(x02, y0,'With FBK','fontsize',fs)
text(x03, y0,'No FBK','fontsize',fs)
text(x04, y0,'With FBK','fontsize',fs)


%First bracket
text(x(1)-0.22,y0-0.003,'__','fontsize',fs-2)
text(x(1)-0.24,y0,'|','fontsize',fs-2)
text(x(1)+0.99,y0-0.003,'__','fontsize',fs-2)
text(x(1)+1.17,y0,'|','fontsize',fs-2)

%Second bracket
text(x(2)-0.22,y0-0.003,'__','fontsize',fs-2)
text(x(2)-0.25,y0,'|','fontsize',fs-2)
text(x(2)+1.19,y0-0.003,'__','fontsize',fs-2)
text(x(2)+1.37,y0,'|','fontsize',fs-2)

%Third bracket
text(x(3)-0.22,y0-0.003,'__','fontsize',fs-2)
text(x(3)-0.24,y0,'|','fontsize',fs-2)
text(x(3)+0.99,y0-0.003,'__','fontsize',fs-2)
text(x(3)+1.17,y0,'|','fontsize',fs-2)

%Fourth bracket
text(x(4)-0.22,y0-0.003,'__','fontsize',fs-2)
text(x(4)-0.25,y0,'|','fontsize',fs-2)
text(x(4)+1.19,y0-0.003,'__','fontsize',fs-2)
text(x(4)+1.37,y0,'|','fontsize',fs-2)

%Attention texts
x01=(xbar(2)+xbar(3))/2-0.4;
x02=(xbar(6)+xbar(7))/2-0.3;
x=[x01 x02];
y0=-0.19;
text(x01-0.2, y0,'NORMAL','fontsize',fs)
text(x02, y0,'HIGH','fontsize',fs)

%Attention brackets
l=0.5; %line length 
a=1.7; %shift factors
b=2.3;
c=0;
for i=1:2
    text(x(i)-l-c,y0+-0.003,'__','fontsize',fs-2)
    text(x(i)-l-0.03-c,y0,'|','fontsize',fs-2)
    text(x(i)+a*l+0.2-c,y0-0.003,'__','fontsize',fs-2)
    text(x(i)+b*l-c+0.09,y0,'|','fontsize',fs-2)
    if mod(i,2)==0
        c=0;
    else
        c=+0.02;
    end
end

if printOn==1
    print(h2,'-depsc','-painters',[path 'figresynchefficiency.eps']);
end

%% FIGURE 4 - Delta MA and delta SD

h3=figure(3);
clf(3);
lwidth=1;
msize=3;
fsize = 9;
fsizelegend = 8;
angle=110; 
density=80;
width=0.6;
fig_size_cm=[20 10];
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters');
set(gcf,'PaperSize',fig_size_cm,'paperposition',[0 0 fig_size_cm]);

% Specific color Setup 
indColors1=[1 1 13 13 27 27 35 35]; 
brillo=brilloInicial;
for i=[1 3 5 7]
    col1=brighten(ametrine(NColors),-brillo);
    col(i,:)=col1(indColors1(i),:);
    col(i+1,:)=col1(indColors1(i),:);
    brillo=brillo-deltaBrillo; 
end

% SUBPLOT - Delta NMA
lineaPunteada=repmat([2.34 -2.34],1,4); %Linear regression (digitalization of fig.2 Repp, panel A 2003) - Rodri 
subplot(11,11,[1:5 12:16 23:27 34:38 45:49 56:60 67:71 78:82 89:93 100:104])

% subplot(1,2,1)
d=0.2;
xbar=[1-d 2-d 3-d 4-d 5+d 6+d 7+d 8+d];
T=[];
M=[];
delta=pos-pre;
T1=squeeze(delta(1,1,:,1,1));
T2=squeeze(delta(1,1,:,1,2));
T3=squeeze(delta(1,1,:,2,1));
T4=squeeze(delta(1,1,:,2,2));
T5=squeeze(delta(1,2,:,1,1));
T6=squeeze(delta(1,2,:,1,2));
T7=squeeze(delta(1,2,:,2,1));
T8=squeeze(delta(1,2,:,2,2));
T9=squeeze(delta(2,1,:,1,1));
T10=squeeze(delta(2,1,:,1,2));
T11=squeeze(delta(2,1,:,2,1));
T12=squeeze(delta(2,1,:,2,2));
T13=squeeze(delta(2,2,:,1,1));
T14=squeeze(delta(2,2,:,1,2));
T15=squeeze(delta(2,2,:,2,1));
T16=squeeze(delta(2,2,:,2,2));
T=[[T1;T9],[T2;T10],[T3;T11],[T4;T12],[T5;T13],[T6;T14],[T7;T15],[T8;T16]];

q=sqrt(NSujetos*2);
for barrai=1:8
    ptp(barrai)=bar(xbar(barrai),mean(T(:,barrai)),'facecolor',col(barrai,:));
    hold on
    h(barrai)=bar(xbar(barrai),lineaPunteada(barrai),'facecolor','none','linestyle','--');
    errorbar(xbar(barrai),mean(T(:,barrai)),std(T(:,barrai))/q,'-','Color',col(barrai,:),'linewidth',lwidth)
    
end

%Scatter plot
for barrai=1:8
    hold on
    for i=1:NSujetos*2
        plot(xbar(barrai)+0.25,T(i,barrai),'o','markerfacecolor','w','markersize',2.5,'linewidth',1,'markeredgecolor','k');
        hhh=plot(xbar(barrai)+0.25,T(i,barrai),'o','markeredgecolor','k','markersize',2.5,'linewidth',1);
        drawnow;
        hhh.MarkerHandle.LineWidth = 0.1;
    end
end

%Hatched bars
hPatch = findobj(h, 'Type', 'bar');
hatchfill2(hPatch,'HatchAngle',angle,'HatchDensity',density,'HatchColor','k','HatchLineWidth',width);

%axis and title
text(-0.145,1.07,'A','units','normalized','fontsize',fsize)
ylim([-60 55]) 
xlim([0 9])
ylabel('{\it \DeltaMA} (ms)','fontsize',fsize,'interpreter','tex')
set(gca,'xtick',[],'xticklabel',[],'fontsize',fsize)
set(gca,'box','off','xcolor','k')

%Significance factors and asteriks
text(0.5,1.07,'attention','fontsize',fsizelegend,'units','normalized')
text(0.5,1.03,'feedback','fontsize',fsizelegend,'units','normalized')
text(0.5,0.99,'pert. sign','fontsize',fsizelegend,'units','normalized')
text(0.7,1.06,'*','fontsize',17,'units','normalized')
text(0.7,1.02,'*','fontsize',17,'units','normalized')
text(0.7,0.98,'***','fontsize',17,'units','normalized')

%vector for ANOVA 4
M=[T(:,1);T(:,2);T(:,3);T(:,4);T(:,5);T(:,6);T(:,7);T(:,8)]; 
sujS=[repmat((1:NSujetos*2)',4,1);repmat((1:NSujetos*2)'+NSujetos*2,4,1)];
atencionA=[repmat(1,1,NSujetos*8)';repmat(2,1,NSujetos*8)'];
feedbackB=[repmat(1,1,NSujetos*4)';repmat(2,1,NSujetos*4)';repmat(1,1,NSujetos*4)';repmat(2,1,NSujetos*4)'];
pertD=repmat([repmat(1,1,NSujetos*2)';repmat(2,1,NSujetos*2)'],4,1);
nesting=[0 0 0 0 ; 0 0 0 0 ; 1 0 0 0 ; 0 0 0 0 ]; 
[p2,table2,stats2,terms2]=anovan(M,{atencionA feedbackB sujS pertD},'random',3,'nested',nesting,'varnames',{'Aten' 'Feed' 'Sujs' 'Pert'},'model','full','display','off');

%Test between bars and expected value
TT=[];
for i=1:8
    if mod(i,2)==1
        ind=1;
    else
        ind=2;
    end
    TT(:,i)=T(:,i)-lineaPunteada(ind);    
end

%Texts p values
% text(0.8,0.95,['Aten p=' num2str(round(1000*table2{2,7})/1000)],'units','normalized','fontsize',13)
% text(0.8,0.9,['Feed p=' num2str(round(1000*table2{3,7})/1000)],'units','normalized','fontsize',13)
% text(0.8,0.85,['Pert p=' num2str(round(1000*table2{5,7})/1000)],'units','normalized','fontsize',13)
% text(0.8,0.8,['Aten*Feed p=' num2str(round(1000*table2{6,7})/1000)],'units','normalized','fontsize',13)
% text(0.8,0.75,['Aten*Pert p=' num2str(round(1000*table2{7,7})/1000)],'units','normalized','fontsize',13)
% text(0.8,0.7,['Feed*Pert p=' num2str(round(1000*table2{9,7})/1000)],'units','normalized','fontsize',13)
% text(0.8,0.65,['Aten*Feed*Pert p=' num2str(round(1000*table2{11,7})/1000)],'units','normalized','fontsize',13)

set(gca,'xtick',xbar,'xticklabel',{'Neg','Pos','Neg','Pos','Neg','Pos','Neg','Pos'},'fontsize',fsize)

%Feedback texts
x01=(xbar(1)+xbar(2))/2-0.65;
x02=(xbar(3)+xbar(4))/2-0.75;
x03=(xbar(5)+xbar(6))/2-0.65;
x04=(xbar(7)+xbar(8))/2-0.75;
x=[x01 x02 x03 x04];
y0=-72;
fs=8;
text(x01, y0,'No FBK','fontsize',fs)
text(x02, y0,'With FBK','fontsize',fs)
text(x03, y0,'No FBK','fontsize',fs)
text(x04, y0,'With FBK','fontsize',fs)



%First bracket
text(x(1)-0.165,y0,'_','fontsize',fs-2)
text(x(1)-0.21,y0+0.03,'|','fontsize',fs-2)
text(x(1)+1.35,y0,'_','fontsize',fs-2)
text(x(1)+1.47,y0+0.03,'|','fontsize',fs-2)

%Second bracket
text(x(2)-0.16,y0,'_','fontsize',fs-2)
text(x(2)-0.195,y0+0.03,'|','fontsize',fs-2)
text(x(2)+1.6,y0,'_','fontsize',fs-2)
text(x(2)+1.72,y0+0.03,'|','fontsize',fs-2)

%Third bracket
text(x(3)-0.165,y0,'_','fontsize',fs-2)
text(x(3)-0.21,y0+0.03,'|','fontsize',fs-2)
text(x(3)+1.35,y0,'_','fontsize',fs-2)
text(x(3)+1.47,y0+0.03,'|','fontsize',fs-2)

%Fourth bracket
text(x(4)-0.16,y0,'_','fontsize',fs-2)
text(x(4)-0.195,y0+0.03,'|','fontsize',fs-2)
text(x(4)+1.6,y0,'_','fontsize',fs-2)
text(x(4)+1.72,y0+0.03,'|','fontsize',fs-2)

%Attention texts
x01=(xbar(2)+xbar(3))/2-0.8;
x02=(xbar(6)+xbar(7))/2-0.5;
x=[x01 x02];
y0=-79;
text(x01, y0,'NORMAL','fontsize',fs)
text(x02, y0,'HIGH','fontsize',fs)

%Attention brackets
d=0.45; %line length 
a=1.7; %shift factors
b=1.95;
c=0;
e=0;
for i=1:2
    text(x(i)-d-c,y0,'__','fontsize',fs-2)
    text(x(i)-d-0.035-c,y0+0.07,'|','fontsize',fs-2)
    text(x(i)+a-c-e,y0,'__','fontsize',fs-2)
    text(x(i)+b-c-e+0.03,y0+0.07,'|','fontsize',fs-2)
    if mod(i,2)==0
        c=0;
        e=0;
    else
        c=+0.05;
        e=0.6;
    end
end

% SUBPLOT - Delta SD 
subplot(11,11,[1:5 12:16 23:27 34:38 45:49 56:60 67:71 78:82 89:93 100:104]+6)
lineaPunteada2=repmat([-1.56 +1.56],1,4); %Linear regression (digitalization of fig.2 Repp, panel B 2003) - Rodri 
d=0.2;
xbar=[1-d 2-d 3-d 4-d 5+d 6+d 7+d 8+d];
T=[];
M=[];

delta=posSTD-preSTD;
T1=squeeze(delta(1,1,:,1,1));
T2=squeeze(delta(1,1,:,1,2));
T3=squeeze(delta(1,1,:,2,1));
T4=squeeze(delta(1,1,:,2,2));
T5=squeeze(delta(1,2,:,1,1));
T6=squeeze(delta(1,2,:,1,2));
T7=squeeze(delta(1,2,:,2,1));
T8=squeeze(delta(1,2,:,2,2));
T9=squeeze(delta(2,1,:,1,1));
T10=squeeze(delta(2,1,:,1,2));
T11=squeeze(delta(2,1,:,2,1));
T12=squeeze(delta(2,1,:,2,2));
T13=squeeze(delta(2,2,:,1,1));
T14=squeeze(delta(2,2,:,1,2));
T15=squeeze(delta(2,2,:,2,1));
T16=squeeze(delta(2,2,:,2,2));
T=[[T1;T9],[T2;T10],[T3;T11],[T4;T12],[T5;T13],[T6;T14],[T7;T15],[T8;T16]];

q=sqrt(NSujetos*2);
for barrai=1:8
    ptp(barrai)=bar(xbar(barrai),mean(T(:,barrai)),'facecolor',col(barrai,:));
    hold on
    h(barrai)=bar(xbar(barrai),lineaPunteada2(barrai),'facecolor','none','linestyle','--');
    errorbar(xbar(barrai),mean(T(:,barrai)),std(T(:,barrai))/q,'-','Color',col(barrai,:),'linewidth',lwidth)
end

%Scatter plot
for barrai=1:8
    hold on
    for i=1:NSujetos*2
        plot(xbar(barrai)+0.25,T(i,barrai),'o','markerfacecolor','w','markersize',2.5,'linewidth',1,'markeredgecolor','k');
        hhh=plot(xbar(barrai)+0.25,T(i,barrai),'o','markeredgecolor','k','markersize',2.5,'linewidth',1);
        drawnow;
        hhh.MarkerHandle.LineWidth = 0.1;
    end
end

%Hatched bars
hPatch = findobj(h, 'Type', 'bar');
hatchfill2(hPatch,'HatchAngle',angle,'HatchDensity',density,'HatchColor','k','HatchLineWidth',width);

%axis and title
text(-0.145,1.07,'B','units','normalized','fontsize',fsize)
ylim([-20 17]) 
xlim([0 9])
ylabel('{\it \DeltaSD} (ms)','fontsize',fsize,'interpreter','tex')
set(gca,'xtick',xbar,'xticklabel',{'Neg','Pos','Neg','Pos','Neg','Pos','Neg','Pos'},'fontsize',fsize)

%Significance factors and asteriks
text(0.5,1.07,'attention   ns','fontsize',fsizelegend,'units','normalized')
text(0.5,1.03,'feedback   ns','fontsize',fsizelegend,'units','normalized')
text(0.5,0.99,'pert. sign','fontsize',fsizelegend,'units','normalized')
text(0.7,0.98,'***','fontsize',17,'units','normalized')

%Feedback texts
x01=(xbar(1)+xbar(2))/2-0.65;
x02=(xbar(3)+xbar(4))/2-0.75;
x03=(xbar(5)+xbar(6))/2-0.65;
x04=(xbar(7)+xbar(8))/2-0.75;
x=[x01 x02 x03 x04];
y0=-24;
fs=8;
text(x01, y0,'No FBK','fontsize',fs)
text(x02, y0,'With FBK','fontsize',fs)
text(x03, y0,'No FBK','fontsize',fs)
text(x04, y0,'With FBK','fontsize',fs)

%First bracket
text(x(1)-0.165,y0,'_','fontsize',fs-2)
text(x(1)-0.21,y0+0.03,'|','fontsize',fs-2)
text(x(1)+1.35,y0,'_','fontsize',fs-2)
text(x(1)+1.47,y0+0.03,'|','fontsize',fs-2)

%Second bracket
text(x(2)-0.16,y0,'_','fontsize',fs-2)
text(x(2)-0.195,y0+0.03,'|','fontsize',fs-2)
text(x(2)+1.6,y0,'_','fontsize',fs-2)
text(x(2)+1.72,y0+0.03,'|','fontsize',fs-2)

%Third bracket
text(x(3)-0.165,y0,'_','fontsize',fs-2)
text(x(3)-0.21,y0+0.03,'|','fontsize',fs-2)
text(x(3)+1.35,y0,'_','fontsize',fs-2)
text(x(3)+1.47,y0+0.03,'|','fontsize',fs-2)

%Fourth bracket
text(x(4)-0.16,y0,'_','fontsize',fs-2)
text(x(4)-0.195,y0+0.03,'|','fontsize',fs-2)
text(x(4)+1.6,y0,'_','fontsize',fs-2)
text(x(4)+1.72,y0+0.03,'|','fontsize',fs-2)

%Attention texts
x01=(xbar(2)+xbar(3))/2-0.8;
x02=(xbar(6)+xbar(7))/2-0.5;
x=[x01 x02];
y0=-26.2;
text(x01, y0,'NORMAL','fontsize',fs)
text(x02, y0,'HIGH','fontsize',fs)

%Attention brackets
d=0.45; %line length 
a=1.7; %shift factors
b=1.95;
c=0;
e=0;
for i=1:2
    text(x(i)-d-c,y0,'__','fontsize',fs-2)
    text(x(i)-d-0.035-c,y0+0.07,'|','fontsize',fs-2)
    text(x(i)+a-c-e,y0,'__','fontsize',fs-2)
    text(x(i)+b-c-e+0.03,y0+0.07,'|','fontsize',fs-2)
    if mod(i,2)==0
        c=0;
        e=0;
    else
        c=+0.05;
        e=0.6;
    end
end
set(gca,'box','off')

%vector for ANOVA 4
M=[T(:,1);T(:,2);T(:,3);T(:,4);T(:,5);T(:,6);T(:,7);T(:,8)]; 
sujS=[repmat((1:NSujetos*2)',4,1);repmat((1:NSujetos*2)'+NSujetos*2,4,1)];
atencionA=[repmat(1,1,NSujetos*8)';repmat(2,1,NSujetos*8)'];
feedbackB=[repmat(1,1,NSujetos*4)';repmat(2,1,NSujetos*4)';repmat(1,1,NSujetos*4)';repmat(2,1,NSujetos*4)'];
pertD=repmat([repmat(1,1,NSujetos*2)';repmat(2,1,NSujetos*2)'],4,1);
nesting=[0 0 0 0 ; 0 0 0 0 ; 1 0 0 0 ; 0 0 0 0 ]; 
[p2,table2,stats2,terms2]=anovan(M,{atencionA feedbackB sujS pertD},'random',3,'nested',nesting,'varnames',{'Aten' 'Feed' 'Sujs' 'Pert'},'model','full','display','off');

%Test between bars and expected value
TT=[];
for i=1:8
    if mod(i,2)==1
        ind=1;
    else
        ind=2;
    end
    TT(:,i)=T(:,i)-lineaPunteada2(ind);    
end

%Print
if printOn==1
    print(h3,'-depsc','-painters',[path 'figDeltaMA.eps']);
end

%%
