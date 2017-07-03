clear;
%%%%%%%%%%%%%%%%%%%%%%%% Getting parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters = [0.4500    1.6674         0   18.6510   32.1000    0.8000    0.0300    0.1200    0.3020    0.0300   60.8500    0.0150 ...
           0.0520   16.0794    0.0800    2.7589    0.4593         0    0.8000];
deltaTF = parameters(12);
constants = CalculateConstantsJY();
FK506=1;
states = CalculateSteadyStatesJY(parameters,FK506);
% constant = struct;
% constant.initialCdt = [0.1,states.switchOff(2),states.switchOff(3),0,0,FK506]; %ori F=0.1
% constant.initial_Switch_Cdt = GetSwitchState(parameters,FK506);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% General Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 2, vary betaRK,deltaTF, keep yNK constant
n = 9;
timespan = [0,2000];
 
switch n        
    case 2 %vary betaRK,deltaTF, keep yNK constant
        healthydeltaTF = parameters(12); %0.05
        deltaTFstart = 0.0001;
        deltaTFincrement = 0.0001;
        deltaTFend = 0.004;
        deltaTFvalue = deltaTFstart:deltaTFincrement:deltaTFend;
        TotalRun = length(deltaTFvalue);
        Result = struct;
        Result.yNK = parameters(3);
        Result.betaRK = [180 300];
        Result.deltaTF = deltaTFvalue;
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.ErrorDeltaTF = [];
        progressbar(0,0);
        for index = 1:length(Result.betaRK)   
            clear parameters(18) %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            Fungus = zeros(1,TotalRun);
            Inflam = zeros(1,TotalRun); 
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for deltaTF = deltaTFstart:deltaTFincrement:deltaTFend
                try
                    clear parameters(12) modelData
                    parameters(12) = deltaTF;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5...
                                    (timespan,[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',deltaTF);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF deltaTF];
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(Result.betaRK),[])
        end
        figure,
        for index = 1:length(Result.betaRK)
            if length(Result.betaRK)>1
                subplot(2,length(Result.betaRK)/2,index)
            end
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            [hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,double(Result.Inflam{index,1}>1e-6));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
            h(3) = line(Result.deltaTF,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');         
            h(3).Color = 'k';
            h(4) = vline(healthydeltaTF,'c','Healthy');
%             h(5) = line(Result.deltaTF,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));  
%             h(5).LineStyle = '--';
           

            for i = 1:size(Result.ErrorDeltaTF)
                vline(Result.ErrorDeltaTF(i),'r');
            end
            
            xlabel('deltaTF','FontSize',18);
            title(['yNK: ' num2str(parameters(3)) ' betaRK: ' num2str(Result.betaRK(index))],'FontSize',16);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus','FontSize',18) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',18) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;  
            grid
        end
        
    case 3 %vary yNK keep betaRK constant
        deltaTFstart = 0.01;
        deltaTFincrement = 0.005;
        deltaTFend = 0.4;
        deltaTFvalue = deltaTFstart:deltaTFincrement:deltaTFend;
        TotalRun = length(deltaTFvalue);
        Result = struct;
        Result.yNK = [4,8,12,16];
        Result.betaRK = parameters(18);
        Result.deltaTF = deltaTFvalue;
        Result.Fungus = cell(length(Result.yNK),1);
        Result.Inflam = cell(length(Result.yNK),1);
        Result.Rswitch = cell(length(Result.yNK),1);

        progressbar(0,0);
        for index = 1:length(Result.yNK)
            clear parameters(3) %parameters(3) = yNK
            parameters(3) = Result.yNK(index);
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for deltaTF = deltaTFstart:deltaTFincrement:deltaTFend
                try
                    clear parameters(12) modelData
                    parameters(12) = deltaTF;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    if modelData(end,5)<1e-6
                        Rswitch(counter) = 0;
                    end
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',deltaTF);
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            mid = 0.5*(max(Fungus)+min(Fungus));
            Result.Rswitch{index,1} = Rswitch.*mid;
            progressbar(index/length(Result.betaRK),[])
        end
        figure,
        for index = 1:length(Result.yNK)
            subplot(2,2,index)
            [hAx,hLine1,hLine2] = plotyy(deltaTFvalue,Result.Fungus{index,1},deltaTFvalue,Result.Inflam{index,1});
            h(3) = line(deltaTFvalue,Result.Rswitch{index,1}, 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            xlabel('deltaTF')
            title(['yNK: ' num2str(Result.yNK(index)) ' betaRK: ' num2str(Result.betaRK)],'FontSize',14);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus','FontSize',13) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',13) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end

    case 4 %vary deltaTF range but keeping yNK and betaRK constant
        deltaTFstart = [0.00095;0.008];
        deltaTFincrement = [0.0002;0.002];
        deltaTFend = [0.01;0.3];
        Result = struct;
        Result.yNK = parameters(3);
        Result.betaRK = parameters(18);
        Result.deltaTF = cell(length(deltaTFstart),1);
        for j = 1:length(deltaTFstart)
            Result.deltaTF{j,1}=deltaTFstart(j):deltaTFincrement(j):deltaTFend(j);
            TotalRun(j) = length(Result.deltaTF{j,1});
        end
        Result.Fungus = cell(length(deltaTFstart),1);
        Result.Inflam = cell(length(deltaTFstart),1);
        Result.Rswitch = cell(length(deltaTFstart),1);
        progressbar(0,0);
        for index = 1:length(deltaTFstart)
            Fungus = ones(1,TotalRun(index));
            Inflam = ones(1,TotalRun(index));
            Rswitch = ones(1,TotalRun(index));
            counter = 1;
            for deltaTF = Result.deltaTF{index,1}
                try
                    clear parameters(12) modelData
                    parameters(12) = deltaTF;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,constant.initialCdt,constant.initial_Switch_Cdt,parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun(index))
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',deltaTF);
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(deltaTFstart),[])
        end
        figure,
        for index = 1:length(deltaTFstart)
            subplot(1,2,index)
            [hAx,hLine1,hLine2] = plotyy(Result.deltaTF{index,1},double(Result.Fungus{index,1}>1e-6),Result.deltaTF{index,1},Result.Inflam{index,1});
            %h(3) = line(yNRvalue,double(Result.Rswitch{index,1}>0), 'Parent', hAx(1), 'Marker','.');
            h(3) = line(Result.deltaTF{index,1},0.1*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');         
            h(3).Color = 'k';
            h(5) = line(Result.deltaTF{index,1},0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));  
            h(5).LineStyle = '--';

            xlabel('deltaTF')
            title(['yNK: ' num2str(Result.yNK), ' betaRK: ' num2str(Result.betaRK)],'FontSize',14);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus','FontSize',13) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',13) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end
        progressbar(1,1)
        
   case 5 %keep deltaTF,yNK constant, vary yNR and betaRK
        yNRstart = 0.00;
        yNRincrement = 0.05;
        yNRend = 1;
        yNRvalue = yNRstart:yNRincrement:yNRend;
        TotalRun = length(yNRvalue);
        Result = struct;
        Result.yNK = parameters(3);
        Result.betaRK = [0.05 1 8 20];
        Result.deltaTF = parameters(12);
        Result.yNR = yNRvalue;
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);

        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for yNR = yNRstart:yNRincrement:yNRend
                try
                    clear parameters(7) modelData
                    parameters(7) = yNR;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
%                     if modelData(end,5)<1e-6
%                         Rswitch(counter) = 0;
%                     end
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when yNR = %d',yNR);
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
%             mid = 0.5;
%             Result.Rswitch{index,1} = Rswitch.*mid;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(Result.betaRK),[])
        end
        figure,
        healthyyNR = 0.005;
        for index = 1:length(Result.betaRK)
            subplot(2,length(Result.betaRK)/2,index)
            [hAx,hLine1,hLine2] = plotyy(yNRvalue,double(Result.Fungus{index,1}>1e-6),yNRvalue,Result.Inflam{index,1});
            %h(3) = line(yNRvalue,double(Result.Rswitch{index,1}>0), 'Parent', hAx(1), 'Marker','.');
            h(3) = line(yNRvalue,0.1*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');         
            h(3).Color = 'k';
            h(4) = vline(healthyyNR,'c','Healthy');
            h(5) = line(yNRvalue,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));  
            h(5).LineStyle = '--';
            xlabel('yNR')
            title(['yNK: ' num2str(Result.yNK) ' betaRK: ' num2str(Result.betaRK(index)), ' deltaTF: ' num2str(Result.deltaTF)],'FontSize',14);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus','FontSize',13) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',13) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end
       
        
       case 6 %keep deltaTF,yNK constant, vary kRN and betaRK
        kRNstart = 5;
        kRNincrement = 0.5;
        kRNend = 45;
        kRNvalue = kRNstart:kRNincrement:kRNend;
        TotalRun = length(kRNvalue);
        Result = struct;
        Result.yNK = parameters(3);
        Result.betaRK = [10 30 100 180];
        Result.deltaTF = parameters(12);
        Result.kRN = kRNvalue;
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);

        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for kRN = kRNstart:kRNincrement:kRNend
                try
                    clear parameters(5) modelData
                    parameters(5) = kRN;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
%                     if modelData(end,5)<1e-6
%                         Rswitch(counter) = 0;
%                     end
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when kRN = %d',kRN);
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
%             mid = 0.5;
%             Result.Rswitch{index,1} = Rswitch.*mid;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(Result.betaRK),[])
        end
        figure,
        healthykRN= 35.3561;
        for index = 1:length(Result.betaRK)
            subplot(2,length(Result.betaRK)/2,index)
            [hAx,hLine1,hLine2] = plotyy(kRNvalue,double(Result.Fungus{index,1}>1e-6),kRNvalue,Result.Inflam{index,1});
            %h(3) = line(yNRvalue,double(Result.Rswitch{index,1}>0), 'Parent', hAx(1), 'Marker','.');
            h(3) = line(kRNvalue,0.1*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');         
            h(3).Color = 'k';
            h(4) = vline(healthykRN,'c','Normal');
            h(5) = line(kRNvalue,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));  
            h(5).LineStyle = '--';
            xlabel('kRN')
            title(['yNK: ' num2str(Result.yNK) ' betaRK: ' num2str(Result.betaRK(index)), ' deltaTF: ' num2str(Result.deltaTF)],'FontSize',14);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus','FontSize',13) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',13) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end
           
        
      case 7 %for a constant betaRK and different yNK, how does that look
        yNKstart = 0;
        yNKincrement = 2;
        yNKend = 80;
        yNKvalue = yNKstart:yNKincrement:yNKend;
        TotalRun = length(yNKvalue);
        Result = struct;
        Result.yNK = yNKvalue;
        Result.betaRK = [0.001 0.01 0.1 1];
        Result.deltaTF = parameters(12);
        Result.kRN = parameters(5);
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);

        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for yNK = yNKstart:yNKincrement:yNKend
                try
                    clear parameters(3) modelData
                    parameters(3) = yNK;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when kRN = %d',kRN);
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(Result.betaRK),[])
        end
        figure,
        for index = 1:length(Result.betaRK)
            subplot(2,length(Result.betaRK)/2,index)
            [hAx,hLine1,hLine2] = plotyy(Result.yNK,double(Result.Fungus{index,1}>1e-6),Result.yNK,double(Result.Inflam{index,1}>1e-6));
            h(3) = line(Result.yNK,0.1*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');         
            h(3).Color = 'k';
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
%             h(5) = line(yNKvalue,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));  
%             h(5).LineStyle = '--';
            xlabel('yNK','FontSize',16)
            title([' betaRK: ' num2str(Result.betaRK(index)), ' deltaTF: ' num2str(Result.deltaTF)],'FontSize',16);
            ylabel('Fungus','FontSize',16);
            ylabel(hAx(1),'Fungus','FontSize',16) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',16) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end

            case 8 %constant parameters, different initial condition
        Fstart = 0;
        Fincrement = 0.4;
        Fend = 10;
        Fvalue = Fstart:Fincrement:Fend;
        TotalRun = length(Fvalue);
        Result = struct;
        Result.initialF = Fvalue;
        Result.betaRK = [5 10 120 200];
        Result.deltaTF = parameters(12);
        Result.kRN = parameters(5);
        Result.yNK = parameters(3);
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Init_switch = GetSwitchState(parameters,FK506);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for initialF = Fstart:Fincrement:Fend
                try
                    clear modelData 
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,[initialF,states.switchOff(2),states.switchOff(3),0,0,FK506],Init_switch,parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when initialF = %d',initialF);
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(Result.betaRK),[])
        end
        figure,
        for index = 1:length(Result.betaRK)
            subplot(2,length(Result.betaRK)/2,index)
            [hAx,hLine1,hLine2] = plotyy(Result.initialF,double(Result.Fungus{index,1}>1e-6),Result.initialF,double(Result.Inflam{index,1}>1e-4));
            h(3) = line(Result.initialF,0.1*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');         
            h(3).Color = 'k';
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
%             h(5) = line(yNKvalue,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));  
%             h(5).LineStyle = '--';
            xlabel('Initial Fungus Value','FontSize',16)
            title([' betaRK: ' num2str(Result.betaRK(index)), ' deltaTF: ' num2str(Result.deltaTF)],'FontSize',16);
            ylabel('Fungus','FontSize',16);
            ylabel(hAx(1),'Fungus','FontSize',16) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',16) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end
        
    case 9 %vary yNK and betaRK, keep deltaTF constant
        betaRKstart = 0.01;
        betaRKincrement = 1;
        betaRKend = 12;
        betaRKvalue = betaRKstart:betaRKincrement:betaRKend;
        TotalRun = length(betaRKvalue);
        Result = struct;
        Result.yNK = [0.05,0.2,0.6,0.8];
        Result.betaRK = betaRKvalue;
        Result.deltaTF = parameters(12);
        Result.Fungus = cell(length(Result.yNK),1);
        Result.Inflam = cell(length(Result.yNK),1);
        Result.Rswitch = cell(length(Result.yNK),1);

        progressbar(0,0);
        for index = 1:length(Result.yNK)
            clear parameters(3) %parameters(3) = yNK
            parameters(3) = Result.yNK(index);
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for betaRK = betaRKstart:betaRKincrement:betaRKend
                try
                    clear parameters(18) modelData
                    parameters(18) = betaRK;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    if modelData(end,5)<1e-6
                        Rswitch(counter) = 0;
                    end
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',deltaTF);
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            mid = 0.5*(max(Fungus)+min(Fungus));
            Result.Rswitch{index,1} = Rswitch.*mid;
            progressbar(index/length(Result.yNK),[])
        end
        figure,
        for index = 1:length(Result.yNK)
            subplot(2,2,index)
            [hAx,hLine1,hLine2] = plotyy(Result.betaRK,double(Result.Fungus{index,1}>1e-6),Result.betaRK,double(Result.Inflam{index,1}>1e-4));
            h(3) = line(Result.betaRK,0.1*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');         
            h(3).Color = 'k';
            xlabel('betaRK')
            title(['yNK: ' num2str(Result.yNK(index)) ' deltaTF: ' num2str(Result.deltaTF)],'FontSize',14);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus','FontSize',13) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',13) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end

        
end
