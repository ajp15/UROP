%Define setting and data for optimisation 
function constants = CalculateConstantsJY()
 
%%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%
constants.ignoreFirstNPar = 0;
constants.annealingSteps = 50;%number of annealing steps for the SA optimistations
constants.optSteps= 50;%number of moves per annealing step for the SA optimisation
%constants.ICs = [3,.1,1,.1,0,0]; %default initial conditions
constants.ignorePar = [1,3,6,10,13,15,18];%parameters to ignore from the optimisations
constants.numRanParameters = 1;%number of iterations for the MC optimisation
constants.error=0.1;%error for costfunction data points
 
constants.lw=8;%line width
constants.fs=28;%font size
constants.ms=15;%marker size
%%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%
 
%fungal growth data 
constants.data.Fgrowth.t = [0,1.5,3,4.5];
fungusLevel = [5, 2.18, 1.36, .452; 5, 2.97, 1.83, .517];
constants.data.Fgrowth.F_MOI3_withM = mean(fungusLevel);

%TNFalpha data
constants.data.TNFalpha.t = [0, 1, 2, 3, 4, 5, 6];
Untreated = [30.63, 29.93, 57.90, 91, 159.29, 184.77, 151.45];
Control = [22.36, 17.75, 12.55, 29.97, 48.82, 49.56, 47.31];
constants.data.TNFalpha.UT = Untreated-Control;
constants.data.TNFalpha.UT = constants.data.TNFalpha.UT - constants.data.TNFalpha.UT(1);

%RCAN1 expression change
constants.data.RCAN1.t = [0, 1, 2, 3, 4, 5, 6];
constants.data.RCAN1.RCAN1levels = [6.0, 5.0, 12.5, 17.5, 27.5, 45.0, 50.0]; %to 1e4
Ratio_RCAN1ToTotalRCAN1 = [0.19, 0.15, 0.17, 0.19, 0.18, 0.20, 0.18];
constants.data.RCAN1.TotalRCAN1levels = constants.data.RCAN1.RCAN1levels./Ratio_RCAN1ToTotalRCAN1;
constants.data.RCAN1.pRCAN1levels = [0.45, 0.425, 0.39, 0.275, 0.25, 0.16, 0.20].*constants.data.RCAN1.TotalRCAN1levels;
constants.data.RCAN1.pRCAN1levels = constants.data.RCAN1.pRCAN1levels - constants.data.RCAN1.pRCAN1levels(1);

%activated NFAT data (initial MOI = 5, #cells = 1e5 )
constants.data.actNFAT.t=[0,0.5,1,1.5,2];
constants.data.actNFAT.UT=[0.082050682, 0.381740934, 0.583135379, 0.705501924,0.885125655];

%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%
 
 
 
end
