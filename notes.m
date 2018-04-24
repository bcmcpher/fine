%% new features / fxns to add to toolbox
% Brent McPherson
%

%% TODO:
%
% streamline render fxns
% - finish / check some functionality
% - better handle individual / multiple inputs in edge render
% - set defualt view to minmax(x, y, z) coords + 10%
%
% make sure volume fxns can take a loaded .nii or a string
%
% feTractProfilePairedConnections - ensure i -> j orientation of profile
%
% fnTractProfileTensor - check functionality
% - add upper diagonal x node output of profiles
%
% fnTractProfileModules - check optimization
%
% fnVirtualLesionModules - failed on 101431; see if it can be fixed
%
% fnPlotModulePoints - plot individuals by color or group average w/ error bars
%
% create classified fiber structure Dan wants
%
% parallel pool fxn
% - see if this is actually necessary - maybe a local problem?
%
% change all display to sprintf
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
% fix ALL documentation
%
% link network functions
% - add other metrics?
% - create separate fxn for computing individual edge similarities?
%
% built in matlab hierarchical clustering algorithms / plots
% linkage / cluster / dendrogram / clusterdata / cophenet / inconsistent
%
% fxn estimating neural conduction speed / bitrate for all edges
%
% summary fxn(s?) of Bassett's network tools
%
% simple shape analysis of cortical node centers; geometric morphometric analysis
% - procrustes alignment
%
% fxn to virtual lesion streamlines w/ both terminations in 1 ROI
%
% better default axes in virtual lesion measures plot 
%
% feCreatePairedConnections 
% - simplify how fibers / weights are stored? too big a change? used later?
%
% finish feCreatePairedConnectionsFromSurface
%
% figure out if single fxn for making matrix field is workable
%
% feTractProfilePairedConnections
% - check if the connections are too short for a reasonable profile (?)
%
% compute profile from ENCODE directly? 
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

% sounds interesting, but not useful
dtiExtrapolateVolumeFromCrossectionArea

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

% compare profiles?
dtiComputeCurveDiff

% apparently from SPM
mrAnatComputeMutualInfo

%% bvalue normalization

% load HCP data
bval = dlmread('~/mrtrix3/dwi.bvals');

% thresholds
bvn = 100;
bv0 = 50;

% grab unique and indices
[ bvals.unique, ~, bvals.uindex ] = unique(bval);

% if bval is <= bv0, assume it's a b-zero
bvals.unique(bvals.unique <= bv0) = 0;

% round to nearest bvn
bvals.unique = round(bvals.unique ./ bvn) * bvn;
bvals.valnorm = bvals.unique( bvals.uindex );
%dlmwrite('dwi.bvals', bvals.valnorm, 'delimiter', ' ');
 
