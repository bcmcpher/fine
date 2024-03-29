%% new features / fxns to add to toolbox
% Brent McPherson
%

%% TODO:
%
% print # of nodes / streamlines to user earlier in fnCreateEdges
%
% export to json graph
% - link network edges?
% - null models?
% - statistics?
%
% better check / parse names if passed during creation
% - will probably handle in specific app fxns - a bl problem
%
% store network stats in object
% - just the binary stats?
% - I doubt the names can stay fixed
% 
% figure out how to store community assignments
% - array of node length with optional labels?
% - optionally check for stored assignments when building modules / plot
%
% add profile tube to render function
%
% built in matlab hierarchical clustering algorithms / plots
% linkage / cluster / dendrogram / clusterdata / cophenet / inconsistent
% - move over from cca_aging (that may be MDS fxn, not sure if that's useful here...)
%
% append to link networks?
% - I'm not sure this is worth the effort
%
% import / store classification.mat
% - basic crossection info on what edges make a fascicle
% - plot of network with edges contributing to fascicles highlighted over
%
% create measures on diagonal self-connections
%
% merge fnPlotModulePoints and fnPlotModuleGroups - they're nearly the same
%
% change all display to sprintf (?)
%
% check all the plot / volume fxn once the majority of above are done
% - they have mostly minor changes / tweaks
% - comments / documentation / examples
%
% make sure volume fxns can take a loaded .nii or a string
%
% fix ALL function names - for consistency / clarity
%
% fix ALL documentation - for consistency / clarity
%
% fxn estimating neural conduction speed / bitrate for all edges
%
% summary fxn(s?) of Bassett's network tools
%
% dual color diagonal adjacency plot - mean / sd in different colors
%
% better default axes in virtual lesion measures plot 
%
% create a way to use surfaces as parc - new assignment works well for it
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
 
