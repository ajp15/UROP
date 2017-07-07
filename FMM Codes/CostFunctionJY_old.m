%Calculates cost, the degree to which the model results deviate from the experimental results
%least squared difference
%the 'cost' is what is minimised by the optimisation
%if plotData==1, then the results of the fits are plotted
function cost = CostFunctionJY(model,parameters,constants,plotData)
 
    plotData = 1;
    parametersInput = parameters;
 
    cost = 0;
 
    if nargin < 4
        %default to not plotting the fits
        plotData=0;
    end
 
    %function that calculates the cost given the model value and the
    %experimental value at each time point
    function cost = CalculateCostJY(constants,times,modelValue,experimentalValue)
        dataError       = experimentalValue*0.1;% 10% marigin of error
        least           = (modelValue - experimentalValue)./dataError;
        cost            = least.*least/length(times); %make sure that more data points do not bias results
    end
 
    %function that matches the model results to the experimental results 
    %and adds the cost of that fit to the total cost
    function cost = AddDataToCost(cost,times,data,modelData,indexToMatch)
        i = 1;
        for(t = times)
            %fuzzy match to find closest time index
            a=abs(modelData(:,1)-t);
            timeIndex = find(a==min(a));
 
            modelValue     = modelData(timeIndex,indexToMatch);%2nd index is F level
            experimentalValue    = data(i);
            cost = cost + CalculateCostJY(constants,times,modelValue,experimentalValue);
 
            i=i+1;
        end
    end
    %to calculate initial condition
    parameters = parametersInput;
    a_temp = (parameters(4)*parameters(7)*parameters(5))/(parameters(2)*parameters(6));
    b = 1+ (parameters(4)/parameters(2));
    c = - 1;
    NFs0   = ( -b + sqrt( b^2 - (4*a_temp*c)))*(1/(2*a_temp));
    RCAN10 = NFs0*parameters(5)/parameters(6);
    %constants.ICs = [3e5,NFs0,RCAN10,0,0,0]; %for model calculation
    IC_FungalData = [5e6,NFs0,RCAN10,0,0,0];
    ICs = [3e5,NFs0,RCAN10,0,0,0];
    initial_Switch_Cdt=GetSwitchState(parameters,0);   
    [modelData_F,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,6],IC_FungalData,initial_Switch_Cdt,parameters);
    [modelData_general,~] = RecursiveHybrid_EventLocationMethod_Fungusv5([0,6],ICs,initial_Switch_Cdt,parameters);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Many different fits corresponding to different sets of experimental data %%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
%    %%%%%%%%%%%%%%%%%%%%%%
%    % Fungal growth data %
%    %%%%%%%%%%%%%%%%%%%%%%
% 
    %%%%
    %with M
    %%%%
    times = constants.data.Fgrowth.t;%F times data
    data = constants.data.Fgrowth.F_MOI3_withM;%F values at these times
    
    %Fit to fungal growth data
    cost = AddDataToCost(cost,times,data,modelData_F,2);
 
    if(plotData==1)%we want to plat the fits
        figure();
        plot(modelData_F(:,1),modelData_F(:,2), 'b', times, data, 'bo','MarkerSize', constants.ms,'LineWidth',constants.lw)
        xlabel('Time[h]', 'FontSize', constants.fs);
        ylabel('F MOI','FontSize', constants.fs);
        title('Fungal Growth data','FontSize', constants.fs)
        hold on;
        %print('-dpdf',strcat(model,'_fungalGrowth_',num2str(CostFunctionJY(model,parameters,constants),'%0.0f')));
    end
 
 
 
    %%%%%%%%%%%%%%%%%
    % TNFalpha data %
    %%%%%%%%%%%%%%%%%
 
    %%%%
    %untreated(w/o FK506)
    %%%%
    %parameters = parametersInput;
    times = constants.data.TNFalpha.t;%T times data
    data = constants.data.TNFalpha.UT;%T values at these times
 
    %Fit to T growth data
    cost = AddDataToCost(cost,times,data,modelData_general,6);
 
    if(plotData==1)
        states = CalculateSteadyStatesJY(parameters)
        figure();
        plot(modelData_general(:,1),modelData_general(:,6), 'b', times, data, 'bo','MarkerSize', constants.ms,'LineWidth',constants.lw)
        hold on;
        plot([0,6],[states.switchOff(5),states.switchOff(5)], ':b', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2)
        plot([0,6],[states.switchOn(5),states.switchOn(5)], ':g', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2)
        xlabel('Time[h]', 'FontSize', constants.fs);
        ylabel('TNFalpha','FontSize', constants.fs);
        title('TNFalpha time course','FontSize', constants.fs)
      %hold on;
    end
 
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % non phosphorylated RCAN1 data %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    %%%%
    %untreated
    %%%%
    times = constants.data.RCAN1.t;%RCAN1 times data
    data = constants.data.RCAN1.RCAN1levels;%RCAN1 fold-change values at these times
    %data=data*RCAN10;%as we have fold-change data (that means RCAN10 here is the 
                     %*actual value* but exp data that we get is fold
                     %change from control. therefore, to get the actual
                     %data, we need to multiply fold change(data) by RCAN10
                 
     %Fit to RCAN1 growth data
    i = 1;
    for(t = times)
        %fuzzy match to find closest time index
        a=abs(modelData_F(:,1)-t);
        timeIndex = find(a==min(a));
 
        modelValue     = modelData_F(timeIndex,4)+modelData_F(timeIndex,5);
        %we calculate this because the experimental data we had is actually
        %total RCAN1 level (ie sum of R and Rp)
        experimentalValue    = data(i);
        cost = cost + CalculateCostJY(constants,times,modelValue,experimentalValue);
        %calculate accummulative cost
        i=i+1;
    end
 
    if(plotData==1)
        states = CalculateSteadyStatesJY(parameters)
 
        figure();
        plot(modelData_F(:,1),modelData_F(:,4)+modelData_F(:,5), 'b', times, data, 'bo','MarkerSize', constants.ms,'LineWidth',constants.lw)
        hold on
        plot(modelData_F(:,1),modelData_F(:,4), 'c',modelData_F(:,1),modelData_F(:,5), 'r' ,'MarkerSize', constants.ms,'LineWidth',constants.lw)
        hold on;
        hline(states.switchOff(3),'c','R SS off');
        hline(states.switchOn(3),'c','R SS n');
        hline(states.switchOff(4),'r','Rp SS off');
        hline(states.switchOn(4),'r','Rp SS n');
        hold on;
        plot([0,7],[parameters(14),parameters(14)], ':k', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2)
        xlabel('Time[h]', 'FontSize', constants.fs);
        ylabel('RCAN1 total','FontSize', constants.fs);
        title('RCAN1 time course','FontSize', constants.fs)
%         print('-dpdf',strcat(model,'_RCAN1data_',num2str(CostFunctionJY(model,parameters,constants),'%0.0f')));
    end
    
    %%%%%%%%%%%%
    % NFs data %
    %%%%%%%%%%%%
 
    %%%%%%%%%%%
    %untreated%
    %%%%%%%%%%%
 
 
    parameters = parametersInput;
    times = constants.data.actNFAT.t;%initial time levels of NFs in the nucleus
    data = constants.data.actNFAT.UT;% 
    cost = AddDataToCost(cost,times,data,modelData_general,3);
  
    if(plotData==1)
        states = CalculateSteadyStatesJY(parameters);
        figure;
        plot(modelData_general(:,1),modelData_general(:,3), 'b', times, data, 'bo','MarkerSize', constants.ms,'LineWidth',constants.lw)
        hold on;
        plot([0,6],[states.switchOff(2),states.switchOff(2)], ':b', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2)
        plot([0,6],[states.switchOn(2),states.switchOn(2)], ':g', 'MarkerSize', constants.ms,'LineWidth',constants.lw/2)
        xlabel('Time[h]', 'FontSize', constants.fs);
        ylabel('NFs','FontSize', constants.fs);
        title('NFs time course','FontSize', constants.fs)
%         print('-dpdf',strcat(model,'_NFs_',num2str(CostFunctionJY(model,parameters,constants),'%0.0f')));
        %hold on; 
    end
 
 
    
 
end
