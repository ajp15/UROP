clear;
%%%%%%%%%%%%%%%%%%%%%%%% Getting parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters = GetMayParameter;
deltaTFwoFK = parameters(12); %0.05
yNRwoFK = parameters(7); %0.0010
kRNwoFK = parameters(5); %35.3561
deltaRwoFK = parameters(6);
deltaRpwoFK = parameters(10);
betaNRwoFK = parameters(9);
kappaRwoFK = parameters(8);
kappaNwoFK = parameters(2);
constants = CalculateConstantsJY();
FK506=1;
states = CalculateSteadyStatesJY(parameters,FK506);
% constant = struct;
% constant.firstinitialCdt = [0.1,states.switchOff(2),states.switchOff(3),0,0,FK506]; %ori F=0.1
% constant.firstinitial_Switch_Cdt = GetSwitchState(parameters,FK506);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% General Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 2, vary betaRK,deltaTF, keep yNK constant
n = 11;

switch n
    case 2 %vary betaRK,deltaTF, keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [120];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        deltaTFstart = 0.03;
        deltaTFincrement = 0.05;
        deltaTFend = 1;
        newdeltaTFvalue = deltaTFstart:deltaTFincrement:deltaTFend;
        Result.deltaTF = newdeltaTFvalue;
        TotalRun = length(newdeltaTFvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(12) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(12) = deltaTFwoFK;
            %parameters(12) = 0.0007; %neutropenia
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,36],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newdeltaTF = deltaTFstart:deltaTFincrement:deltaTFend
                try
                    clear parameters(12) modelData
                    parameters(12) = newdeltaTF;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',newyNR);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newyNR];
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
            h(4) = vline(deltaTFwoFK,'c','No Treatment');
            %            h(5) = vline(0.0007,'r','No Treatment'); %neutropenia
            %             h(5) = line(Result.deltaTF,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));
            %             h(5).LineStyle = '--';
            
            xlabel('deltaTF after treatment','FontSize',18);
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
        
    case 3 %vary betaRK,yNR,keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [5,20,120,250];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        Result.deltaTF = deltaTFwoFK;
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        newyNRstart = 0.00005;
        newyNRincrement = 0.0005;
        newyNRend = 0.02;
        newyNRvalue = newyNRstart:newyNRincrement:newyNRend;
        TotalRun = length(newyNRvalue);
        Result.yNR = newyNRvalue;
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(7) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(7) = yNRwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,720],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newyNR = newyNRstart:newyNRincrement:newyNRend
                try
                    clear parameters(7) modelData
                    parameters(7) = newyNR;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,1000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',newyNR);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newyNR];
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(Result.betaRK),[])
        end
        
        figure,
        counter = 1;
        for index = [1,3,4]%length(Result.betaRK)
%             if length(Result.betaRK)>1
%                 subplot(2,length(Result.betaRK)/2,index)
%             end
            subplot(1,3,counter);
            counter=counter+1;
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            
            [hAx,hLine1,hLine2] = plotyy(Result.yNR,double(Result.Fungus{index,1}>1e-6),Result.yNR,double(Result.Inflam{index,1}>1e-6));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1]) %h(3) = line(yNRvalue,double(Result.Rswitch{index,1}>0), 'Parent', hAx(1), 'Marker','.');
            h(3) = line(Result.yNR,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            h(4) = vline(yNRwoFK,'r','No Treatment');
            %             h(5) = line(Result.yNR,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));
            %             h(5).LineStyle = '--';
            xlabel('yNR after treatment','FontSize',18)
            title(['yNK: ' num2str(Result.yNK) ' betaRK: ' num2str(Result.betaRK(index)), ' deltaTF: ' num2str(Result.deltaTF)],'FontSize',16);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus (F)','FontSize',18) % left y-axis
            ylabel(hAx(2),'TNFalpha (T)','FontSize',18) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end
        
    case 4 %vary betaRK,kRN,keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [5,120,200];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        Result.deltaTF = deltaTFwoFK;
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        newkRNstart = 0;
        newkRNincrement = 1;
        newkRNend = 40;
        newkRNvalue = newkRNstart:newkRNincrement:newkRNend;
        TotalRun = length(newkRNvalue);
        Result.kRN = newkRNvalue;
        TreatmentTime = 72;
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(5) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(5) = kRNwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,TreatmentTime],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newkRN = newkRNstart:newkRNincrement:newkRNend
                try
                    clear parameters(5) modelData
                    parameters(5) = newkRN;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when kRN = %d',newkRN);
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
            if length(Result.betaRK)==3
                subplot(1,3,index)
            elseif length(Result.betaRK)>1
                subplot(2,length(Result.betaRK)/2,index)
            end
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            [hAx,hLine1,hLine2] = plotyy(Result.kRN,double(Result.Fungus{index,1}>5e-0),Result.kRN,double(Result.Inflam{index,1}>5e-0));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1.4]) %h(3) = line(yNRvalue,double(Result.Rswitch{index,1}>0), 'Parent', hAx(1), 'Marker','.');
            h(3) = line(Result.kRN,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            h(4) = vline(kRNwoFK,'r','No Treatment');
            %             h(5) = line(Result.kRN,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));
            %             h(5).LineStyle = '--';
            xlabel('kRN after treatment','FontSize',18)
            title(['yNK: ' num2str(Result.yNK) ' betaRK: ' num2str(Result.betaRK(index)), ' deltaTF: ' num2str(Result.deltaTF)],'FontSize',20);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus (F)','FontSize',18) % left y-axis
            ylabel(hAx(2),'TNFalpha (T)','FontSize',18) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid;
        end
        
    case 5 %vary betaRK,deltaTF, keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [0.1 120];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        deltaTFstart = 0.0001;
        deltaTFincrement = 0.0002;
        deltaTFend = 0.01;
        newdeltaTFvalue = deltaTFstart:deltaTFincrement:deltaTFend;
        Result.deltaTF = newdeltaTFvalue;
        TotalRun = length(newdeltaTFvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(12) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(12) = deltaTFwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,72],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newdeltaTF = deltaTFstart:deltaTFincrement:deltaTFend
                try
                    clear parameters(12) modelData
                    parameters(12) = newdeltaTF;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,928],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',newyNR);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newyNR];
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
            [hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,double(Result.Inflam{index,1}>1e-0));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
            h(3) = line(Result.deltaTF,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            %            h(4) = vline(deltaTFwoFK,'c','No Treatment');
            h(5) = vline(0.0007,'r','No Treatment'); %neutropenia
            %             h(5) = line(Result.deltaTF,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));
            %             h(5).LineStyle = '--';
            
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
        
    case 10 %vary betaRK,deltaTF, keep yNK constant
        %define the combination of betaRK and yNK
        load('temp.mat')
        Fungus_short = Result.Fungus;
        Inflam_short = Result.Inflam;
        Rswitch = Result.Rswitch;
        
        Result = struct;
        Result.betaRK = [0.1 120];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        deltaTFstart = 0.0001;
        deltaTFincrement = 0.0002;
        deltaTFend = 0.01;
        newdeltaTFvalue = deltaTFstart:deltaTFincrement:deltaTFend;
        Result.deltaTF = newdeltaTFvalue;
        TotalRun = length(newdeltaTFvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(12) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            %parameters(12) = deltaTFwoFK;
            parameters(12) = 0.0007;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,1440],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newdeltaTF = deltaTFstart:deltaTFincrement:deltaTFend
                try
                    clear parameters(12) modelData
                    parameters(12) = newdeltaTF;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaTF = %d',newyNR);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newyNR];
                    continue
                end
            end
            Result.Fungus{index,1} = Fungus;
            Result.Inflam{index,1} = Inflam;
            Result.Rswitch{index,1} = Rswitch;
            progressbar(index/length(Result.betaRK),[])
        end
        
        figure,
        for index = 1:length(Result.betaRK)+1
            
            subplot(3,1,index)
            if index == 1
                [hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Fungus_short{1,1}>1e-6),Result.deltaTF,double(Inflam_short{1,1}>1e-0));
                set(hAx(1),'YLim',[0 1.4])
                set(hAx(2),'YLim',[-1 1])
                h(3) = line(Result.deltaTF,0.5*double(Rswitch>1e-6), 'Parent', hAx(1), 'Marker','.');
                h(3).Color = 'k';
                %            h(4) = vline(deltaTFwoFK,'c','No Treatment');
                h(5) = vline(0.0007,'r','No Treatment'); %neutropenia
                %             h(5) = line(Result.deltaTF,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));
                %             h(5).LineStyle = '--';
                
                xlabel('deltaTF after treatment at 36th hour','FontSize',15);
                title('yNK: 4.5 betaRK:0.1','FontSize',16);
                ylabel('Fungus');
                ylabel(hAx(1),'Fungus','FontSize',18) % left y-axis
                ylabel(hAx(2),'TNFalpha','FontSize',18) % right y-axis
                hLine1.Marker = '*';
                hLine2.Marker = 'x';
                hLine1.LineWidth = 3;
                hLine2.LineWidth = 3;
                grid
            else
                
                %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
                [hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index-1,1}>1e-6),Result.deltaTF,double(Result.Inflam{index-1,1}>1e-0));
                set(hAx(1),'YLim',[0 1.4])
                set(hAx(2),'YLim',[-1 1])
                h(3) = line(Result.deltaTF,0.5*double(Result.Rswitch{index-1,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
                h(3).Color = 'k';
                %            h(4) = vline(deltaTFwoFK,'c','No Treatment');
                h(5) = vline(0.0007,'r','No Treatment'); %neutropenia
                %             h(5) = line(Result.deltaTF,0.5.*max(Result.Inflam{index,1}).*double(Result.Inflam{index,1}>1e-6), 'Parent', hAx(2));
                %             h(5).LineStyle = '--';
                
                xlabel('deltaTF after treatment at 1440th hour','FontSize',15);
                title(['yNK: ' num2str(parameters(3)) ' betaRK: ' num2str(Result.betaRK(index-1))],'FontSize',16);
                ylabel('Fungus');
                ylabel(hAx(1),'Fungus','FontSize',18) % left y-axis
                ylabel(hAx(2),'TNFalpha','FontSize',18) % right y-axis
                hLine1.Marker = '*';
                hLine2.Marker = 'x';
                hLine1.LineWidth = 3;
                hLine2.LineWidth = 3;
                grid
            end
        end
        
    case 6 %vary betaRK,deltaR, keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [5 120 200];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        deltaRstart = 0.02;
        deltaRincrement = 0.02;
        deltaRend = 0.5;
        newdeltaRvalue = deltaRstart:deltaRincrement:deltaRend;
        Result.deltaR = newdeltaRvalue;
        TotalRun = length(newdeltaRvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(6) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(6) = deltaRwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,360],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newdeltaR = deltaRstart:deltaRincrement:deltaRend
                try
                    clear parameters(6) modelData
                    parameters(6) = newdeltaR;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaR = %d',newdeltaR);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newdeltaR];
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
            
            if length(Result.betaRK)==3
                subplot(1,3,index);
            elseif length(Result.betaRK)>1
                subplot(2,length(Result.betaRK)/2,index)
            end
            
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            [hAx,hLine1,hLine2] = plotyy(Result.deltaR,double(Result.Fungus{index,1}>1e-6),Result.deltaR,double(Result.Inflam{index,1}>1e-6));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
            h(3) = line(Result.deltaR,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            h(5) = vline(deltaRwoFK,'r','No Treatment'); %neutropenia
            
            xlabel('deltaR after treatment','FontSize',18);
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
        
    case 7 %vary betaRK,deltaR, keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [5 120 200];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        deltaRpstart = 0.06;
        deltaRpincrement = 0.02;
        deltaRpend = 0.8;
        newdeltaRpvalue = deltaRpstart:deltaRpincrement:deltaRpend;
        Result.deltaRp = newdeltaRpvalue;
        TotalRun = length(newdeltaRpvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(10) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(10) = deltaRpwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,360],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newdeltaRp = deltaRpstart:deltaRpincrement:deltaRpend
                try
                    clear parameters(10) modelData
                    parameters(10) = newdeltaRp;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaRp = %d',newdeltaRp);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newdeltaRp];
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
            
            if length(Result.betaRK)==3
                subplot(1,3,index);
            elseif length(Result.betaRK)>1
                subplot(2,length(Result.betaRK)/2,index)
            end
            
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            [hAx,hLine1,hLine2] = plotyy(Result.deltaRp,double(Result.Fungus{index,1}>1e-0),Result.deltaRp,double(Result.Inflam{index,1}>5e1));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
            h(3) = line(Result.deltaRp,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            h(5) = vline(deltaRpwoFK,'r','No Treatment'); %neutropenia
            
            xlabel('deltaRp after treatment','FontSize',18);
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
        
         case 8 %vary betaRK,deltaR, keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [5 120 200];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        betaNRstart = 0.1;
        betaNRincrement = 0.1;
        betaNRend = 1.8;
        newbetaNRvalue = betaNRstart:betaNRincrement:betaNRend;
        Result.betaNR = newbetaNRvalue;
        TotalRun = length(newbetaNRvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(9) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(9) = betaNRwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,360],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newbetaNR = betaNRstart:betaNRincrement:betaNRend
                try
                    clear parameters(9) modelData
                    parameters(9) = newbetaNR;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when deltaR = %d',newdeltaR);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newdeltaR];
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
            
            if length(Result.betaRK)==3
                subplot(1,3,index);
            elseif length(Result.betaRK)>1
                subplot(2,length(Result.betaRK)/2,index)
            end
            
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            [hAx,hLine1,hLine2] = plotyy(Result.betaNR,double(Result.Fungus{index,1}>1e-6),Result.betaNR,double(Result.Inflam{index,1}>1e-0));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
            h(3) = line(Result.betaNR,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            h(5) = vline(betaNRwoFK,'r','No Treatment'); %neutropenia
            
            xlabel('betaNR after treatment','FontSize',18);
            title(['yNK: ' num2str(parameters(3)) ' betaRK: ' num2str(Result.betaRK(index))],'FontSize',20);
            ylabel('Fungus');
            ylabel(hAx(1),'Fungus','FontSize',18) % left y-axis
            ylabel(hAx(2),'TNFalpha','FontSize',18) % right y-axis
            hLine1.Marker = '*';
            hLine2.Marker = 'x';
            hLine1.LineWidth = 3;
            hLine2.LineWidth = 3;
            grid
        end

          case 9 %vary betaRK,kappaR, keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [5 120 200];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        kappaRstart = 0.0001; %0.001 normal
        kappaRincrement = 0.00005;
        kappaRend = 0.0014;
        newkappaRvalue = kappaRstart:kappaRincrement:kappaRend;
        Result.kappaR = newkappaRvalue;
        TotalRun = length(newkappaRvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(8) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(8) = kappaRwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,360],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newkappaR = kappaRstart:kappaRincrement:kappaRend
                try
                    clear parameters(10) modelData
                    parameters(8) = newkappaR;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when kappaR = %d',newkappaR);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newkappaR];
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
            
            if length(Result.betaRK)==3
                subplot(1,3,index);
            elseif length(Result.betaRK)>1
                subplot(2,length(Result.betaRK)/2,index)
            end
            
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            [hAx,hLine1,hLine2] = plotyy(Result.kappaR,double(Result.Fungus{index,1}>5e-0),Result.kappaR,double(Result.Inflam{index,1}>5e1));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
            h(3) = line(Result.kappaR,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            h(5) = vline(kappaRwoFK,'r','No Treatment'); %neutropenia
            
            xlabel('kappaR after treatment','FontSize',18);
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

         case 11 %vary betaRK,kappaR, keep yNK constant
        %define the combination of betaRK and yNK
        Result = struct;
        Result.betaRK = [5 120 200];
        Result.Fungus = cell(length(Result.betaRK),1);
        Result.Inflam = cell(length(Result.betaRK),1);
        Result.Rswitch = cell(length(Result.betaRK),1);
        Result.yNK = parameters(3);
        IC_data = zeros(length(Result.betaRK),6);
        IC_switch = zeros(length(Result.betaRK),2);
        kappaNstart = 0.1; %1.3161 normal
        kappaNincrement = 0.04;
        kappaNend = 1.4;
        newkappaNvalue = kappaNstart:kappaNincrement:kappaNend;
        Result.kappaN = newkappaNvalue;
        TotalRun = length(newkappaNvalue);
        progressbar(0,0);
        for index = 1:length(Result.betaRK)
            clear parameters(18) parameters(8) BeforeTreatmentData BeforeTreatmentSwitch_state %parameters(18) = betaRK
            parameters(18) = Result.betaRK(index);
            parameters(2) = kappaNwoFK;
            [BeforeTreatmentData,BeforeTreatmentSwitch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,360],[0.1,states.switchOff(2),states.switchOff(3),0,0,FK506],GetSwitchState(parameters,FK506),parameters);
            %BeforeTreatmentSwitch_state = 2xN
            IC_data(index,:) =  BeforeTreatmentData(end,2:7);
            IC_switch(index,:) = transpose(BeforeTreatmentSwitch_state(:,end));
            Fungus = ones(1,TotalRun);
            Inflam = ones(1,TotalRun);
            Rswitch = ones(1,TotalRun);
            counter = 1;
            for newkappaN = kappaNstart:kappaNincrement:kappaNend
                try
                    clear parameters(10) modelData
                    parameters(2) = newkappaN;
                    [modelData,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,2000],IC_data(index,:),IC_switch(index,:),parameters);
                    Fungus(counter) = modelData(end,2);
                    Inflam(counter) = modelData(end,6);
                    Rswitch(counter) = modelData(end,5);
                    progressbar([],counter/TotalRun)
                    counter = counter+1;
                catch
                    fprintf('Error when kappaR = %d',newkappaN);
                    Result.ErrorDeltaTF = [Result.ErrorDeltaTF newkappaN];
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
            
            if length(Result.betaRK)==3
                subplot(1,3,index);
            elseif length(Result.betaRK)>1
                subplot(2,length(Result.betaRK)/2,index)
            end
            
            %[hAx,hLine1,hLine2] = plotyy(Result.deltaTF,double(Result.Fungus{index,1}>1e-6),Result.deltaTF,Result.Inflam{index,1});
            [hAx,hLine1,hLine2] = plotyy(Result.kappaN,double(Result.Fungus{index,1}>1e-0),Result.kappaN,double(Result.Inflam{index,1}>1e-0));
            set(hAx(1),'YLim',[0 1.4])
            set(hAx(2),'YLim',[-1 1])
            h(3) = line(Result.kappaN,0.5*double(Result.Rswitch{index,1}>1e-6), 'Parent', hAx(1), 'Marker','.');
            h(3).Color = 'k';
            h(5) = vline(kappaNwoFK,'r','No Treatment'); %neutropenia
            
            xlabel('kappaN after treatment','FontSize',18);
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

        
end

%%


