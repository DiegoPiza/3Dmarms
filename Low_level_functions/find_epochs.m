function epochs = find_epochs( cc, threshold, sign,gap )
% *WAVE*
%
% FIND EPOCHS      locate epochs in a logical array, based on some simple
%                       heuristic parameters
%
% INPUT
% cc - correlation array (1,t)
% threshold - thresholding value, to produce a logical array (sc)
% Edits by DBP 2022
% threshold sign - positive = true (more than >), negative = false (less
% than < )
% gap (in number of samples) allows gaps of gap many samples (positive integer)
% OUTPUT
% epochs - structure with epoch start/end times
%

assert( isvector(cc) == 1, 'vector input required, cc' )
assert( isscalar(threshold), 'scalar input required, threshold' )

% init
epochs = struct([]);

% threshold the logical array
if sign
    L = ( cc > threshold );
elseif ~sign
    L = ( cc <= threshold );
end

d = [true, diff(~L') ~= 0, true];  % TRUE if values change
n = diff(find(d));               % Number of repetitions
Y = repelem(n, n);
Y(~L'==0)=0; %LENGTH OF THE GAPS


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

            if ( epoch_found && Y(ii)>gap )
                
                epochs(epoch_number).end_time = ii;
                epoch_found = 0;
                
            elseif ( epoch_found && Y(ii)<gap )
                
                % do nothing -- allowing gaps of gap sample
                
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
