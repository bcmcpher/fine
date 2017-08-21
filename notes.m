%% new measures to add to edge-wise nodes for more thorough anatomy
% Brent McPherson
%

%% TODO:

% feTractProfilePairedConnections
% - pass dt6 instead of msvol to create multiple profiles
% - check if a labels exists and don't overwrite / add clobber
% - pass nnodes / minnum as fxn level arguments
% - check if the connections are too short for a reasonable profile (?)
%
% feVirtualLesionPairedConnections
% - add 'norm' option for computing virtual lesions
% - check that the results I mean() are supposed to be meaned
%
% feCreatePairedConnections 
% - catch endpoints as ROI data
%
% feCleanPairedConnections
% - provide fxn level access to cleaning parameters

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
% USE CODE HERE TO COMPUTE WEIGHTED DISTANCE
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


