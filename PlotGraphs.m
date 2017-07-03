function PlotGraphs(ModelData,Switch_state,~,parameters,constants,timeWindow,setPrint,plotSwitch,TreatmentTime)
   
    if nargin < 9
        TreatmentTime = 0;
    end
    
    timeStart=timeWindow(1);
    timeEnd=timeWindow(2);

    states = CalculateSteadyStatesJY(parameters,ModelData(1,7));
    
    %for plotting of the switch statuses
    if plotSwitch==1
        Switch_state(1,:)=parameters(14).*1.05.*Switch_state(1,:);
        Switch_state(2,:)=parameters(19).*1.05.*Switch_state(2,:);
    end

    % plot fungal growth data
    figure();
    %semilogy(ModelData(:,1),ModelData(:,3), 'c',ModelData(:,1),ModelData(:,4), 'r',ModelData(:,1),ModelData(:,5), 'k',ModelData(:,1),ModelData(:,6), 'm',ModelData(:,1),ModelData(:,2),'b','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    semilogy(ModelData(:,1),ModelData(:,3), 'c',ModelData(:,1),ModelData(:,4), 'r',ModelData(:,1),ModelData(:,5), 'k','MarkerSize', constants.ms,'LineWidth',3);hold on,
    semilogy(ModelData(:,1),ModelData(:,6), 'm',ModelData(:,1),ModelData(:,2),'b','MarkerSize', constants.ms,'LineWidth',5);
  
    hold on
    %NFs ss
    plot(timeWindow,[parameters(19),parameters(19)], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    %plot(timeWindow,[states.Ron_Non_FK506(2),states.Ron_Non_FK506(2)],'-.c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    %if ModelData(1,7) == 0
        plot(timeWindow,[states.switchOff(2),states.switchOff(2)], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
        plot(timeWindow,[states.switchOn(2),states.switchOn(2)], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     else
%         plot(timeWindow,[states.Ron_Non_FK506(2),states.Ron_Non_FK506(2)], '-.c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%         plot(timeWindow,[states.Roff_Noff_FK506(2),states.Roff_Noff_FK506(2)], '-.c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%         
   % end
    %     plot(timeWindow,[states.Roff_Non_Fnotclc(2),states.Roff_Non_Fnotclc(2)],':b', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    
    %RCan1 ss
    plot(timeWindow,[parameters(14),parameters(14)], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    %if  ModelData(1,7) == 0
        plot(timeWindow,[states.switchOff(3),states.switchOff(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
        plot(timeWindow,[states.switchOn(3),states.switchOn(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     else
%         plot(timeWindow,[states.Ron_Non_FK506(3),states.Ron_Non_FK506(3)], '-.r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%         plot(timeWindow,[states.Roff_Noff_FK506(3),states.Roff_Noff_FK506(3)], '-.r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     end
    
    if plotSwitch==1
        plot(ModelData(:,1),Switch_state(1,:),'y','MarkerSize', constants.ms,'LineWidth',constants.lw/4)
        plot(ModelData(:,1),Switch_state(2,:),'y','MarkerSize', constants.ms,'LineWidth',constants.lw/4)
    end
    
    if TreatmentTime ~= 0
        vline(TreatmentTime,'g:')
    end
    xlabel('Time[h]', 'FontSize', constants.fs);
    ylabel('variables','FontSize', constants.fs);
    axis([timeStart timeEnd 1e-5 2e4]) %max(ModelData(:,6))
    
%     
%     if  ModelData(1,7) == 1
%         title (['deltaTF = ' num2str(parameters(12)) ' yNK = ' num2str(parameters(3)) ' betaRK = ' num2str(parameters(18))],'FontSize',18)
%     end
%     if setPrint==1
%         print('-dpdf',model);
%     end
%     
%     
%     
%     %Phase diagram
%     figure();
%     plot(ModelData(:,3),ModelData(:,4), 'b','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
%     hold on
%     %NFs ss
%     plot([states.switchOff(2),states.switchOff(2)],[min(ModelData(:,4)),max(ModelData(:,4))], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     plot([states.switchOn(2),states.switchOn(2)],[min(ModelData(:,4)),max(ModelData(:,4))], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     plot([parameters(19),parameters(19)],[min(ModelData(:,4)),max(ModelData(:,4))], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     %RCan1 ss
%     plot([min(ModelData(:,3)),max(ModelData(:,3))],[states.switchOff(3),states.switchOff(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     plot([min(ModelData(:,3)),max(ModelData(:,3))],[states.switchOn(3),states.switchOn(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     plot([min(ModelData(:,3)),max(ModelData(:,3))],[parameters(14),parameters(14)], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
%     xlabel('active NFAT (N*)', 'FontSize', constants.fs);
%     ylabel('Total RCAN1','FontSize', constants.fs);
%     axis([min(ModelData(:,3)),max(ModelData(:,3)),min(ModelData(:,4)),max(ModelData(:,4))])
%     
%     if setPrint==1
%         print('-dpdf',strcat('FungalMacrophageModel','phase'));
%     end
%     
%     %Print T against F
%     figure();
%     plot(ModelData(:,2),ModelData(:,6),'b','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
%     xlabel('Fungus', 'FontSize', constants.fs);
%     ylabel('TNFalpha','FontSize', constants.fs);
%     axis([min(ModelData(:,2)),max(ModelData(:,2)),min(ModelData(:,6)),max(ModelData(:,6))])
% 
% 


end

