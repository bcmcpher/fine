%% new measures to add to edge-wise nodes for more thorough anatomy
% Brent McPherson
%

%% TODO:

% 
% store curv/tors parameters in profiles, not a separate field
%
% link network functions
% - add option to skip other similarity measures
% - add average angle between tracts
% - deal with norms between profiles or in a separate fxn?
% - documentation
%
% parallel pool fxn
% - see if this is actually necessary - maybe a local problem?
%
% streamline render fxns
% - finish / check some functionality
% - better handle individual / multiple inputs in edge render
% - set defualt view to minmax(x, y, z) coords + 10%
%
% figure out if single fxn for making matrix field is workable
%
% fill in profiles as nan x nnodes so size can be reliably found 
% -improve fnTractProfileTensor
%
% plot_endpoint_heatmap
% - add option to smooth
% - add option to reslice
%
% fnEstimateLouvainCommunity
% - accept a range of gamma
% - assign default parameters
%
% feTractProfilePairedConnections
% - check if the connections are too short for a reasonable profile (?)
%
% fnBinaryNetworkStats
% - assign / otherwise mark a default value for stats
%
% fnRentianScaling
% - deal with rois structure or make a separate helper to pull that?
% - set up a defualt number of iterations
%
% feCreatePairedConnections 
% - superfiber w/ usable default parameters
% - catch endpoints as ROI data
% - catch standard deviation of length
% - return object that has total assigned counts / descriptives
% - simplify how fibers / weights are stored? too big a change? used later?
% - must find a minimum number of streamlines here (?)
%
% feCleanPairedConnections
% - provide fxn level access to cleaning parameters
% - explicitly pass a minimum number of streamlines for a connection to exist
%
% plotVirtualLesion (LiFE plots)
% - fix overlap line
% - separate handles?
% - better internal notes
%
% fsInflateDestriuex
% - actually make it work
%
% fsInflateDK
% - add option for subcortical / cerebellar?
%
% scripts and workflows are essentially the same - merge somehow
% - drop project specific fxns, but don't loose them...
%
% check all the plot / volume fxn once the majority of above are done
% - they have mostly minor changes / tweaks
% - comments / documentation / examples
%
% fix ALL function names
%
% fxn estimating neural conduction speed / bitrate for all edges
%

%% figure out network shape stats for edges

% load a dt6 to see what some of this stuff is
dt6 = dtiLoadDt6('101431/dt6_xflip/dti48trilin/dt6.mat');

% the actaul easiest way
% Westinshapes - the shape data I want to add
eigVal = dtiGetAllEigValsFromFibers(dt6, fg_class(1));
[ cl, cp, cs ] = dtiComputeWestinShapes(eigVal, 'l1');
% every node of every streamline computed here

% should this be done on a super streamline instead?
% a tract profile that is the WS - see if that's how it's done

% all require dt6
[ eigVec, eigVal ] = dtiEig(dt6);
dtArr = dtiEigComp(eigVec, eigVal);

% relevant but not useful
dtiEigenvaluesFromWestinShapes();

% See what other shape summaries I can pull from something like
% USE CODE HERE TO COMPUTE WEIGHTED DISTANCE FOR A FG PROFILE
[ curv, tors ] = AFQ_ParameterizeTractShape(fg);

% mba code to compute values by streamline
% compute on super-fiber for edge-wise stats?
% see AFQ_ParameterizeTractShape for applying weights from super-fiber
mbaFiberProperties

% compute torsion of single fiber
% average or compute super fiber first?
dtiFiberTorsion

dtiComputeDiffusionPropertiesAlognFG
% has many extra outputs with an included dt6
% I will not rely on dt6 and shouldn't need them for edge shapes, will work around
% but it's a good reference to see what it uses
dtiFiberGroupPropertyWeightedAverage % is internal

% In AFQ is an external tool
frenet2
% computes curvature and torsion for an FG, among other things
% better / other tool to use?

% compute curvature stats from vistasoft on fg group (Hiromasa)
dtiComputeFiberCurvature
dtiComputeFiberCurvatureDistribution
dtiComputeFiberRadiusofCurvatureDistribution

% create super-fiber for some other edge things
dtiComputeSuperFiberRepresentation
dtiComputeSuperFiberMeans
dtiFiberResample
dtiMergeSuperFiberGroups
dtiReorientSuperfiberGroups
dtiReorientFibers
dtiReorderFibersByDistance
dtiResampleFiberGroup

% maybe this can be passed pconn?
dtiFiberProperties
dtiFiberSummary

%% others / extras

% add optional smoothing / output space to endpoint ROI map

% create probability map in MNI space for each edge across subjects, similar to:
AFQ_MakeFGProbabilityMap();

% compute dt6 anat to MNI xform - can save down for group edge probability maps
dtiComputeFgToMNItransform(dt6);

% just remember to check what it does
dtiCheckMotion

% interesting map of tensor properties I've never seen reported
dtiCoherenceMap

% what? cool...
dtiNeuralConductionSpeed
dtiNeuronBitRate


