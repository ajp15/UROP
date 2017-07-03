clear;
%Defining parameters and settings
constants = CalculateConstantsJY(); 
%parameters = GetParameter;
parameters = [0.4500   38.6674         0  386.6740    8.1000    0.2000    2.8000...
    0.0002   80.0000    0.0300   30.8500    0.0500    0.0520...
12.0794    0.0008    2.7589    0.4593    0.9500    0.4000];
ratio = parameters(18)/parameters(3);

timespan = [0,1000];
states = CalculateSteadyStatesJY(parameters,0); %states = CalculateSteadyStatesJY(parameters,FK506);
initialCdt = [1e6,states.switchOff(2),states.switchOff(3),0,0,0];
initial_Switch_Cdt = GetSwitchState(parameters,0);
%parameters(14) = states.switchOn(3)*1.15; 
%parameters(19) = states.switchOn(2)*1.2;
%Integration
%[modelData,Switch_state] = RecursiveHybrid_EventLocationMethod_Fungusv2(timespan,initialCdt,initial_Switch_Cdt,parameters);
timeStart = tic;
[modelData,Switch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5(timespan,initialCdt,initial_Switch_Cdt,parameters);
tElapsed = toc(timeStart)

PlotGraphs(modelData,Switch_state,states,parameters,constants,timespan,0,0) 
%CheckCondition(parameters,yNK_end,betaRK_end)
%openfig('Cdt41v3.fig');
