function [ fe, x, y ] = feTestNiter(fe, outfile)
%% observe change in weights over many LiFE iterations
% Brent McPherson
% 20170509
%
% figure; hold on; for ii = 1:size(y, 1); plot(x, y(ii, :)); end
%

%% load data

%/N/dc2/projects/lifebid/franpest/161223_M

display('Loading the fe structure...');

% load the fe structure
load(fe);

% define colors for iterations
%cols = {[ 0 0 0 ], [ .9 0 0 ], [ 0 0 .9 ], [ 0 .9 0 ], [.3 .7 0 ], [.7 .3 0 ]};

% pull the number of fibers
nFibers = feGet(fe, 'nfibers');

% create empty initial weights
w0 = zeros(nFibers, 1);

% create empty histogram y
y = zeros(75, 100);

%figure;
%hold on;

for rep = 1:75
    
    display(['Repeat #: ' num2str(rep)]);
    
    % fit the weights, provide the previous starting weights
    Niter = 200;
    fe = feSet(fe, 'fit', feFitModel(feGet(fe, 'model'), ...
                                     feGet(fe, 'dsigdemeaned'), 'bbnnls', ...
                                     Niter, 'preconditioner', w0));
    
    % pull weights from fe
    w0 = feGet(fe, 'fiberweights');
    
    % begin the plot
    if rep == 1
        
        display(['Storing Repeat #: ' num2str(rep)]);
        
        % pull bin centers on the first iteration
        [ y(rep, :), x ] = hist(log10(w0(w0 > 0)), 100);
        % plot the first histogram
        %plot(x, y(rep, :), 'color', cols{rep}); 
        save(outfile, 'x', 'y', '-v7.3');
        
    else
        
        display(['Storing Repeat #: ' num2str(rep)]);
        
        % assign values
        [ y(rep, :), x ] = hist(log10(w0(w0 > 0)), x);
        % plot the first histogram
        %plot(x, y(rep, :), 'color', cols{rep});
        save(outfile, 'x', 'y', '-append');
        
    end
    
end

display('Saving final output...');

% save all the outputs
save(outfile, 'fe', '-append');

end
