function epochs = find_epochs( cc, threshold )
% *WAVE*
%
% FIND EPOCHS      locate epochs in a logical array, based on some simple
%                       heuristic parameters
%
% INPUT
% cc - correlation array (1,t)
% threshold - thresholding value, to produce a logical array (sc)
%
% OUTPUT
% epochs - structure with epoch start/end times
%

assert( isvector(cc) == 1, 'vector input required, cc' )
assert( isscalar(threshold), 'scalar input required, threshold' )

% init
epochs = struct([]);

% threshold the logical array
L = ( cc > threshold );

% loop through L
epoch_found = 0;
epoch_number = 0;
for ii = 1:length(L)-1
    
    if L(ii)
        
        if ~epoch_found
            
            epoch_found = 1;
            epoch_number = epoch_number + 1;
            epochs(epoch_number).start_time = ii; % working in INDs
            
        elseif epoch_found
            
            % do nothing
            
        end
        
    elseif ~L(ii)
        
        if ( epoch_found && ~L(ii+1) )
            
            epochs(epoch_number).end_time = ii;
            epoch_found = 0;
            
        elseif ( epoch_found && L(ii+1) )
            
            % do nothing -- allowing gaps of one sample
            
        end
        
    end
    
end

% check output
for ii = 1:length(epochs)
    
    if ~isfield( epochs(ii), 'end_time' )
        epochs(ii).end_time = length(L);
    end
    
end

% calculate summary statistics
for ii = 1:length(epochs)
    
    epochs(ii).epoch_length = epochs(ii).end_time - epochs(ii).start_time;
        
end
