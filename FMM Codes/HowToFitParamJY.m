%JYHee 2016
%Run this code to perform optimisation.
%Begin by providing an initial guess of parameters. 
%     a           = parameters(1);%fungal growth rate
%     kN         = parameters(2);%NFAT activation rate
%     yNK         = parameters(3);%FK506 inhibition rate (of NFAT activation)
%     deltaN      = parameters(4);%NFAT deactivation rate
%     kRN          = parameters(5);%RCAN1 production rate
%     deltaR      = parameters(6);%RCAN1 degredation rate
%     yNR         = parameters(7);%RCAN1 inhibition rate (of NFAT activation)
%     kR         = parameters(8);%RCAN1 production rate
%     betaNR         = parameters(9);%Enhanced NFAT activation rate (by RCAN1-p)
%     deltaRP     = parameters(10);%RCAN1-p degredation rate
%     kT          = parameters(11);%TNFalpha production
%     deltaTF          = parameters(12);%Fungal killing rate
%     deltaT      = parameters(13);%TNFalpha degredation rate
%     Rc           = parameters(14);%RCAN1 switch threshold
%     x           = parameters(15);%RCAN1 switch hill coefficient
%     betaNF         = parameters(16);%NFAT activation enhancement by F
%     betaRF         = parameters(17);%RCAN1-p production enhancement by F
%     betaRK         = parameters(18);%RCAN1-p production enhancement by F
%     Nc         = parameters(19);%Threshold activated NFAT causes inflammation

%There are 2 optimisation methods available:
%1. Monte-Carlo Simulated Annealing (MCSA): slower but more likely to find the parameter set you want
%2. Monte-Carlo (MC): faster but will only really fine tune a good guess of initialParameters.
%You can change the settings for the optimisations by editing the CalculateConstantsJY() file. 
%e.g.
%constants.annealingSteps - defines how precise to make the parameters in the MCSA optimisation routine (the bigger the more precise)
%constants.optSteps - defines how much to search for the best parameter set (the bigger, the better fits you?ll get)
%constants.numRanParameters - same as above but for the MC optimisation (number 2)
%constants.ignorePar - defines the index of the parameters that you don?t want to optimise , i.e. the parameters that you know and want to stay the same 

function HowToFitParamJY() 
clear; close all;

initialParameters = [00.4500   38.6674         0  386.6740    8.1000    0.2000 ...
    2.8000    0.0002   80.0000    0.0300   30.8500    0.0500    0.0520 ...
    12.0794    0.0008    2.7589    0.4593    0.9500    0.4000 ...
];
                 
initialParameters(4) = initialParameters(2)/0.1; 
  
%choose your optimisation method
%optimisation = 1;%MSCA
%optimisation = 2;%MC
 optimisation=1;
 
%now run the optimisation
ControlJY(initialParameters,optimisation)
end 
 
%when the optimisation is running the current minimum cost is printed to the console
%if nothing appears the optimisation is having difficulty improving on initialParameters
 
%when the optimisation is finished graphs will appear in the directory that 
%show how well the parameter fitting has gone and the final parameter set is output to a .txt file in csv format,
%the first number being the cost.


