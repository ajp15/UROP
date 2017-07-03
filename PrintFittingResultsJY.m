%function that prints long-term system dynamics
function yOutput=PrintFittingResultsJY(model,initialParameters,constants)

    parameters = initialParameters;
    
    T=1000;

    constants.cost=1;%CostFunction(model,initialParameters,constants);

    states = CalculateSteadyStatesJY(parameters);

    if(states.switchOff(2)>0)   %to check if initial condition is below 0
        constants.ICs = [0.1,states.switchOff(2),states.switchOff(3),0,0,0];%MOI=3
        %constants.ICs=[F,N,R,Rp,T,FK506]
    else
        error 'IC (states.switchOff(2)) below 0'
    end
    
    function [SwitchState] = GetSwitchState(parameters)
        if states.switchOff(3)>parameters(14)
            Rswitch=1;
        else
            Rswitch=0;
        end
        
        if states.switchOff(2)>parameters(19)
            Nswitch=1;
        else
            Nswitch=0;
        end
        
        SwitchState=[Rswitch,Nswitch];
    end
    
    initial_Switch_Cdt=GetSwitchState(parameters);
    
    timespan=[0,T];
   
    [modelData,Switch_state] = RecursiveHybrid_EventLocationMethod_Fungus(timespan,constants.ICs,initial_Switch_Cdt,parameters,parameters(14),parameters(19));
    %[modelData,Switch_state] = RecursiveHybrid_EventLocationMethod_Fungusv1(timespan,constants.ICs,initial_Switch_Cdt,parameters,parameters(14),parameters(19));

    %for plotting of the switch statuses
    Switch_state(1,:)=30.*Switch_state(1,:);
    Switch_state(2,:)=(0.2).*Switch_state(2,:);

    % plot fungal growth data
    figure();
    semilogy(modelData(:,1),modelData(:,2), 'b', modelData(:,1),modelData(:,3), 'c',modelData(:,1),modelData(:,4), 'r',modelData(:,1),modelData(:,5), 'k',modelData(:,1),modelData(:,6), 'm','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    hold on
    %NFs ss
    plot([0,T],[states.switchOff(2),states.switchOff(2)], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([0,T],[states.switchOn(2),states.switchOn(2)], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([0,T],[initialParameters(19),initialParameters(19)], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    %RCan1 ss
    plot([0,T],[states.switchOff(3),states.switchOff(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([0,T],[states.switchOn(3),states.switchOn(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([0,T],[parameters(14),parameters(14)], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot(modelData(:,1),Switch_state(1,:),'.r','MarkerSize', constants.ms,'LineWidth',constants.lw/4)
    plot(modelData(:,1),Switch_state(2,:),'.c','MarkerSize', constants.ms,'LineWidth',constants.lw/4)
    xlabel('Time[h]', 'FontSize', constants.fs);
    ylabel('variables','FontSize', constants.fs);
    axis([0 T 1e-2 1e3])
    print('-dpdf',model);


    %Phase diagram
    figure();
    plot(modelData(:,3),modelData(:,4), 'b','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    hold on
    %NFs ss
    plot([states.switchOff(2),states.switchOff(2)],[min(modelData(:,4)),max(modelData(:,4))], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([states.switchOn(2),states.switchOn(2)],[min(modelData(:,4)),max(modelData(:,4))], '--c', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([initialParameters(19),initialParameters(19)],[min(modelData(:,4)),max(modelData(:,4))], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    %RCan1 ss
    plot([min(modelData(:,3)),max(modelData(:,3))],[states.switchOff(3),states.switchOff(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([min(modelData(:,3)),max(modelData(:,3))],[states.switchOn(3),states.switchOn(3)], '--r', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    plot([min(modelData(:,3)),max(modelData(:,3))],[initialParameters(14),initialParameters(14)], '--k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/3)
    xlabel('NFAT active', 'FontSize', constants.fs);
    ylabel('RCAN1','FontSize', constants.fs);
    axis([min(modelData(:,3)),max(modelData(:,3)),min(modelData(:,4)),max(modelData(:,4))])
    print('-dpdf',strcat(model,'phase'));

    %Print T against F
    figure();
    plot(modelData(:,2),modelData(:,6),'b','MarkerSize', constants.ms,'LineWidth',constants.lw/2)
    xlabel('F', 'FontSize', constants.fs);
    ylabel('T','FontSize', constants.fs);
    axis([min(modelData(:,2)),max(modelData(:,2)),min(modelData(:,6)),max(modelData(:,6))])


end