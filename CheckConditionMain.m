clear;
%parameters = GetParameter;
parameters= [0.4500    1.6674         0   18.6510   32.1000    0.8000    0.0300    0.1200    0.3020    0.0300   60.8500    0.0150 ...
           0.0520   16.0794    0.0800    2.7589    0.4593         0    0.8000];

%healthydeltaTF =  0.05;
yNK_end = 1;
betaRK_end = 8;
n = 4; 

alpha           = parameters(1);%fungal growth rate
kN         = parameters(2);%NFAT activation rate
yNK         = parameters(3);%FK506 inhibition rate (of NFAT activation)
deltaN      = parameters(4);%NFAT deactivation rate
kRN          = parameters(5);%RCAN1 production rate
deltaR      = parameters(6);%RCAN1 degredation rate
yNR         = parameters(7);%RCAN1 inhibition rate (of NFAT activation)
kR         = parameters(8);%RCAN1 production rate
betaNR         = parameters(9);%Enhanced NFAT activation rate (by RCAN1-p)
deltaRP     = parameters(10);%RCAN1-p degredation rate
kT          = parameters(11);%TNFalpha production
deltaTF          = parameters(12);%Fungal killing rate
deltaT      = parameters(13);%TNFalpha degredation rate
Rc           = parameters(14);%RCAN1 switch threshold
x           = parameters(15);%RCAN1 switch hill coefficien08
betaNF         = parameters(16);%NFAT activation enhancement by F
betaRF         = parameters(17);%RCAN1-p production enhancement by F
betaRK         = parameters(18);%RCAN1-p production enhancement by F
Nc         = parameters(19);%Threshold activated NFAT causes inflammation

clear betaRK yNK

%deltaTF = 0.008;
%yNR = healthy_yNR *1.5;
FK506=1;
row=1;
col=1;
inflammation=[];
yNKincrement=0.1;
bRKincrement=0.5;
betaRKvalue=0:bRKincrement:betaRK_end;
yNKvalue=0:yNKincrement:yNK_end;
SumDifference=zeros(length(yNKvalue),length(betaRKvalue));
PIwithF = zeros(length(yNKvalue),length(betaRKvalue));
PIonly = zeros(length(yNKvalue),length(betaRKvalue));
TrueRatio = zeros(size(PIwithF));
PF = zeros(size(PIwithF));
NoPI_index = [];
PI_index = [];
LC_index = [];

switch n
    case 1 %Focus on PI despite fungal clearance
        for yNK=0:yNKincrement:yNK_end
            for betaRK=0:bRKincrement:betaRK_end
                kR2=kR*(1+betaRK*FK506);
                kN2=kN/(1+yNK*FK506);
                a=(betaNR*kR2*kRN)/(deltaRP*(kR2+deltaR))+(deltaN*yNR*kRN)/(kN2*(kR2+deltaR));
                b=1+deltaN/kN2-(betaNR*kR2*kRN)/(deltaRP*(kR2+deltaR));
                c=-1;
                NFs=(1/(2*a))*(-b+sqrt(b^2-4*a*c));

                Rss=kRN*NFs/(kR2+deltaR);

                if (NFs-Nc)>0 && (Rss-Rc)>0
                    PIonly(row,col)=500; %cannot use 100, too light
                    TrueRatio(row,col) = betaRK/yNK;
                end             
                col=col+1;
            end
           col=1;
           row=row+1;
        end
    Type_colormap = brewermap(16,'RdGy');
        
    case 2 %Persistent fungal growth
        for yNK=0:yNKincrement:yNK_end
            for betaRK=0:bRKincrement:betaRK_end               
                kR2=kR*(1+betaRK*FK506);
                kN2=kN/(1+yNK*FK506);
                Tss = alpha/deltaTF;
                NFs = deltaT*Tss/kT;
                a_Cdt3 = (1-NFs)*kN2*betaNF*betaRF*(kR2+betaNR*kR2*kRN*NFs/deltaRP);
                b2 = betaRF*(kR2+betaNR*kR2*kRN*NFs/deltaRP);
                b3 = betaNF*(kR2+deltaR+betaNR*kR2*kRN*NFs/deltaRP);
                b_Cdt3 = (1-NFs)*kN2*(b2+b3)-deltaN*NFs*kR2*betaRF;
                c_Cdt3 = (1-NFs)*kN2*(kR2+deltaR+betaNR*kR2*kRN*NFs/deltaRP)-deltaN*NFs*(kR2+deltaR+yNR*kRN*NFs);
                Fss = (1/(2*a_Cdt3))*(-b_Cdt3+sqrt(b_Cdt3^2-4*a_Cdt3*c_Cdt3));
                Rss = kRN*NFs/(kR2*(1+betaRF*Fss)+deltaR); 
                Rpss = kR2*(1+betaRF*Fss)*Rss/deltaRP;
                
                if (Fss>0)
                   PF(row,col) = 800; 
                   if (Rss>Rc) && (Fss>0) && (NFs>Nc)
                       PIwithF(row,col) = 500;
                       TrueRatio(row,col) = yNK/betaRK;                   
                   end
                end
                col=col+1;
            end
           col=1;
           row=row+1;
        end
    Type_colormap = brewermap(16,'RdGy');
        
    case 3 %check for both PI, PI+PF, no ratio discrimination 
        for yNK=0:yNKincrement:yNK_end
            for betaRK=0:bRKincrement:betaRK_end
                kR2=kR*(1+betaRK*FK506);
                kN2=kN/(1+yNK*FK506);
                a=(betaNR*kR2*kRN)/(deltaRP*(kR2+deltaR))+(deltaN*yNR*kRN)/(kN2*(kR2+deltaR));
                b=1+deltaN/kN2-(betaNR*kR2*kRN)/(deltaRP*(kR2+deltaR));
                c=-1;
                NFs=(1/(2*a))*(-b+sqrt(b^2-4*a*c));


                Rss=kRN*NFs/(kR2+deltaR);
                Rpss_cdt1 = kR2*Rss/deltaRP;

                if (NFs>Nc) && (Rss>Rc) 
                    %PI(row,col)=500;
                    PIonly(row,col)=500;
                end

                Tss = alpha/deltaTF; 
                NFs_Cdt3 = deltaT*Tss/kT;
                a_Cdt3 = (1-NFs_Cdt3)*kN2*betaNF*betaRF*(kR2+betaNR*kR2*kRN*NFs_Cdt3/deltaRP);
                b2 = betaRF*(kR2+betaNR*kR2*kRN*NFs_Cdt3/deltaRP);
                b3 = betaNF*(kR2+deltaR+betaNR*kR2*kRN*NFs_Cdt3/deltaRP);
                b_Cdt3 = (1-NFs_Cdt3)*kN2*(b2+b3)-deltaN*NFs_Cdt3*kR2*betaRF;
                c_Cdt3 = (1-NFs_Cdt3)*kN2*(kR2+deltaR+betaNR*kR2*kRN*NFs_Cdt3/deltaRP)-deltaN*NFs_Cdt3*(kR2+deltaR+yNR*kRN*NFs_Cdt3);
                Fss = (1/(2*a_Cdt3))*(-b_Cdt3+sqrt(b_Cdt3^2-4*a_Cdt3*c_Cdt3));
                Rss_Cdt3 = kRN*NFs_Cdt3/(kR2*(1+betaRF*Fss)+deltaR); 
                Rpss = kR2*(1+betaRF*Fss)*Rss_Cdt3/deltaRP;

                if Fss>0
                    PF(row,col) = 900;
                    if (Rss_Cdt3>Rc) && (Fss>0) && (NFs_Cdt3>Nc)
                        PIwithF(row,col) = 600;
                    end
                end
                col=col+1;
           end
           col=1;
           row=row+1;
        end
   Type_colormap = brewermap(8,'BuPu');
   
    case 4 %check condition by real simulation 
        clear betaRKvalue yNKvalue SumDifference PIwithF PIonly
        timespan = [0,2000];
        states = CalculateSteadyStatesJY(parameters,FK506);
        constant = struct;
        constant.initialCdt = [0.1,states.switchOff(2),states.switchOff(3),0,0,FK506];
        constant.initial_Switch_Cdt = GetSwitchState(parameters,FK506);
        Result = struct;
        Result.yNK = 0:0.05:1;
        Result.betaRK = 0:0.5:10;
        yNKvalue = Result.yNK;
        betaRKvalue = Result.betaRK;
        Result.param = parameters;
        Fungus = zeros(length(Result.yNK),length(Result.betaRK));
        Inflam = zeros(length(Result.yNK),length(Result.betaRK));
        Rswitch = zeros(length(Result.yNK),length(Result.betaRK));
        Outcome = zeros(length(Result.yNK),length(Result.betaRK),4);
        TotalRun = length(Result.yNK)*length(Result.betaRK);
        progressbar(0);
        counter = 1;
        for yNK = Result.yNK   
            clear parameters(3) %parameters(18) = betaRK
            parameters(3) = yNK;
            for betaRK = Result.betaRK
                try
                    clear parameters(18) modelData
                    parameters(18) = betaRK;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,constant.initialCdt,constant.initial_Switch_Cdt,parameters);
                    Fungus(row,col) = modelData(end,2);
                    Inflam(row,col) = modelData(end,6);
                    Rswitch(row,col) = modelData(end,5);
                    progressbar(counter/TotalRun)
                    counter = counter+1;
                    Outcome(row,col,1) = double((Fungus(row,col)<1e-6)&&(Inflam(row,col)>1e-6)&&(Rswitch(row,col)>1e-6));
                    Outcome(row,col,2) = double((Fungus(row,col)>1e-6)&&(Inflam(row,col)>1e-6)&&(Rswitch(row,col)<1e-6));
                    Outcome(row,col,3) = double((Fungus(row,col)>1e-6)&&(Inflam(row,col)>1e-6)&&(Rswitch(row,col)>1e-6));
                    Outcome(row,col,4) = double((Fungus(row,col)<1e-6)&&(Inflam(row,col)>1e-6)&&(Rswitch(row,col)<1e-6)); %
                    col = col+1;
                catch err
                    col = col+1;  
                    fprintf('Error when yNK = %d, betaRK = %d;\n',yNK,betaRK);
                    disp(getReport(err,'extended'));
                    continue
                end
            end
            col = 1;
            row = row+1;
        end
        Result.Fungus = Fungus;
        Result.Inflam = Inflam;
        Result.Rswitch = Rswitch;
        Result.Outcome = Outcome; 
        PIonly = Outcome(:,:,1);
        PIwithF = Outcome(:,:,3);
        PIwithF_Roff = Outcome(:,:,2);
        if nnz(Outcome(:,:,2))>0
            display('Outcome 2 non zero!');
        end
        Type_colormap = brewermap(16,'BuPu');
        
    case 5 %checking for outcome 2
        for yNK=0:yNKincrement:yNK_end
            for betaRK=0:bRKincrement:betaRK_end  
                kN2=kN/(1+yNK*FK506);
                NFs = alpha*deltaT/(kT*deltaTF);
                Rss = kRN*NFs/deltaR;
                Fss = (1/betaNF)*(deltaN*NFs*(1+yNR*Rss)/(kN2*(1-NFs))-1);
                
                if Rss < Rc && NFs > Nc && Fss > 0
                    PIwithF = 300;
                end
                
            end
        end
        Type_colormap = brewermap(16,'RdGy');
end
    
    PIwithF(isinf(PIwithF))=1000; %remove Inf data

    figure
   ShowPIwithF = imagesc(betaRKvalue,yNKvalue,PIwithF); %both PI and PF
   hold on,
%     ShowRatio = imagesc(betaRKvalue,yNKvalue,TrueRatio);
%     set(ShowRatio,'AlphaData', 1);
%     hold on
%     ShowPF = imagesc(betaRKvalue,yNKvalue,PF);
%     set(ShowPF,'AlphaData', 0.2);
%     hold on,
    ShowPIonly = imagesc(betaRKvalue,yNKvalue,PIonly);
    set(ShowPIonly,'AlphaData', 0.2);
    hold on,
%     ShowPIwithF_Roff= imagesc(betaRKvalue,yNKvalue,PIwithF_Roff);
%     set(ShowPIwithF_Roff,'AlphaData', 0.08);
%     hold on,
    colormap(Type_colormap);
    set(gca,'YDir','normal')
    colorbar
    xlabel('betaRK','FontSize',20); ylabel('yNK','FontSize',20);
    title(['deltaTF = ' num2str(deltaTF), ' yNR = ' num2str(yNR)],'FontSize',16);
    hold on,
    PlotDataPt(deltaTF,yNR)
    if isempty(Result)==0
        save('Condition3General.mat','Result');
    end

