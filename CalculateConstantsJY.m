%Want to optimise ...
function constants = CalculateConstantsJY()
 
 
%%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%
constants.ignoreFirstNPar = 0;
constants.annealingSteps = 1;%number of annealing steps for the SA optimistations
constants.optSteps= 1;%number of moves per annealing step for the SA optimisation
%constants.ICs = [3,.1,1,.1,0,0]; %initial conditions
constants.ignorePar = [1,3,6,10,13,15,18];%parameters to ignore from the optimisations
constants.numRanParameters = 1;%number of iterations for the MC optimisation
constants.error=0.1;%error for costfunction data points
 
constants.lw=8;%line width
constants.fs=28;%font size
constants.ms=15;%marker size
%%%%%%%%%%%%%SETTINGS%%%%%%%%%%%%%%%
 
 
%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%
 
%fungal growth data - assuming M values are approximately constant, track MOI
% constants.data.Fgrowth.t = [0,0.5,1,1.5,2,3,4,5,6];
% constants.data.Fgrowth.F_MOI3_withM = [1,1.368766,1.45487,1.604514,1.516628,3.407959,5.950716,6.306417,15.51664]*3;%with M, MOI of 3
constants.data.Fgrowth.t = [0,1.5,3,4.5];
fungusLevel = [5e6, 2.18e6, 1.36e6, 4.52e5; 5e6, 2.97e6, 1.83e6, 5.17e5];
constants.data.Fgrowth.F_MOI3_withM = mean(fungusLevel);

%TNFalpha data
% constants.data.TNFalpha.t = [0,0.5,1,1.5,2,3,4,5,6];
% constants.data.TNFalpha.UT = [6.526628,3.347712,8.040946,6.051882,7.267049,9.124006,40.78624,144.7928,90.04564]-6.526628;
% constants.data.TNFalpha.FK506 = [0.777828,0.787457,1.60187,4.21619,4.34573,12.57963,30.84656,42.03423,48.1425];

constants.data.TNFalpha.t = [0, 1, 2, 3, 4, 5, 6];
Untreated = [30.63, 29.93, 57.90, 91, 159.29, 184.77, 151.45];
Control = [22.36, 17.75, 12.55, 29.97, 48.82, 49.56, 47.31];
constants.data.TNFalpha.UT = Untreated-Control;
constants.data.TNFalpha.UT = constants.data.TNFalpha.UT - constants.data.TNFalpha.UT(1);

%RCAN1 expression change
%constants.data.RCAN1.t = [0,0.5,1,1.5,2,3,4,5,6];
%constants.data.RCAN1.RCAN1levels = [1,0.845913,0.897191,0.854498, 1.884128, 2.08796, 3.542307,2.878536,3.047512]; %fold change, must be multiplied by SS_R_low
constants.data.RCAN1.t = [0, 1, 2, 3, 4, 5, 6];
constants.data.RCAN1.RCAN1levels = [6.0, 5.0, 12.5, 17.5, 27.5, 45.0, 50.0]; %to 1e4
Ratio_RCAN1ToTotalRCAN1 = [0.19, 0.15, 0.17, 0.19, 0.18, 0.20, 0.18];
constants.data.RCAN1.TotalRCAN1levels = constants.data.RCAN1.RCAN1levels./Ratio_RCAN1ToTotalRCAN1;
constants.data.RCAN1.pRCAN1levels = [0.45, 0.425, 0.39, 0.275, 0.25, 0.16, 0.20].*constants.data.RCAN1.TotalRCAN1levels;
constants.data.RCAN1.pRCAN1levels = constants.data.RCAN1.pRCAN1levels - constants.data.RCAN1.pRCAN1levels(1);


% %activated NFAT data (initial MOI = 5, #cells = 1e5 )
constants.data.actNFAT.t=[0,0.5,1,1.5,2];
constants.data.actNFAT.UT=[0.082050682, 0.381740934, 0.583135379, 0.705501924,0.885125655];
 
% 
% constants.data.actNFAT.t=[0,0.25,0.5,1,1.5,2,3,4,5,6];
% constants.data.actNFAT.UT=[0.358,0.537,0.466,0.447,0.416,0.447,0.381,0.339,0.425,0.368];
% %constants.data.actNFAT.FK506=[0.358,];
% 
%constants.data.actNFAT.t=[0,0.5,1,1.5,2];
%constants.data.actNFAT.UT=[1,1.348,1.567,0.724,0.229];
%constants.data.actNFAT.FK506=[0.066,0.154,0.282,0.193,0.353];
 
%%%%%%%%%%%%%DATA%%%%%%%%%%%%%%%
 
 
 
 
 
end
