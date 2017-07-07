%JYHee 2016
%For manual parameter fitting. Generate the model result without doing
%optimisation and plot the experimental result on the same graph as the
%simulated result
clear;

parameters = [00.4500   38.6674         0  386.6740    8.1000    0.2000 ...
    2.8000    0.0002   80.0000    0.0300   30.8500    0.0500    0.0520 ...
    12.0794    0.0008    2.7589    0.4593    0.9500    0.2000];

%Setting%
model='FungalMacrophageModel';
constants = CalculateConstantsJY(); 
timeWindow=[0,6];
states = CalculateSteadyStatesJY(parameters,0);

%Initial Condition and integration
initialCdt = [5e6,states.switchOff(2),6, 0, 0, 0]; %for model calculation
initial_Switch_Cdt=GetSwitchState(parameters,0);  
[modelData_F,Switch_state] = RecursiveHybrid_EventLocationMethod_Fungusv5(timeWindow,initialCdt,initial_Switch_Cdt,parameters);
initialCdt = [5e6,states.switchOff(2),states.switchOff(3),0,0,0]; 
[modelData_NFs,Switch_stateNFs] = RecursiveHybrid_EventLocationMethod_Fungusv5(timeWindow,initialCdt,initial_Switch_Cdt,parameters);

PrintResults(modelData_F,modelData_NFs,states,parameters,Switch_state,Switch_stateNFs,1)