%Neville J. Boon 13/11/2015
%this document describes how to fit parameters to data for the fungal macrophage model
%it can be run in matlab - just copy paste for an example! (but make sure you are in the directory where all the m-files are)
    
%define your initial guess of the parameters. For example:
%initialParameters = [0.3375,0.5,0,0.005,20,0.5,90.225834,1.522147,0.262957,0.230000,355.185101,0.000000,2.285000,0.674727,5.000000,6.034617,0.175982,0.368105,0.2];
function HowToFitParamJY() 
clear; close all;

initialParameters = [00.4500   38.6674         0  386.6740    8.1000    0.2000 ...
    2.8000    0.0002   80.0000    0.0300   30.8500    0.0500    0.0520 ...
    12.0794    0.0008    2.7589    0.4593    0.9500    0.4000 ...
];
                 
initialParameters(4) = initialParameters(2)/0.1;
%choose what data you want to optimise to
%see examples in the CalculateConstants() file
%there must be a corresponding element in the CostFunction() file to define how to edit the model to the data and what plot will demonstrate visually how well the optimisation has gone
 
%Now choose the settings for the optimisations - edit the CalculateConstants() file. 
%e.g.
%constants.annealingSteps - defines how precise to make the parameters in the MCSA optimisation routine (the bigger the more precise)
%constants.optSteps - defines how much to search for the best parameter set (the bigger, the better fits you?ll get)
%constants.numRanParameters - same as above but for the MC optimisation (number 2)
%constants.ignorePar - defines the index of the parameters that you don?t want to optimise , i.e. the parameters that you know and want to stay the same 
 
%choose your optimisation method
%optimisation = 1;%MCSA: slower but more likely to find the parameter set you want
%optimisation = 2;%MC: faster but will only really fine tune a good guess of initialParameters.
 
%   Columns 1 through 13
% 
%     0.4500   38.6674         0  386.6740    8.1000    0.2000    2.8000    0.0002   80.0000    0.0300   30.8500    0.0350    0.0520
% 
%   Columns 14 through 19
% 
%    12.0794    0.0080    2.7589    0.4593    0.9500    0.2000

optimisation=2;
 
%now run the optimisation
ControlJY(initialParameters,optimisation)
end 
 
%when the optimisation is running the current minimum cost is printed to the console
%if nothing appears the optimisation is having difficulty improving on initialParameters
 
%when the optimisation is finished graphs will appear in the directory that 
%show how well the parameter fitting has gone and the final parameter set is output to a .txt file in csv format, the first number being the cost.


