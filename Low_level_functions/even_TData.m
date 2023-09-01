function TrackingData = even_TData(TrackingData)
    % Ensures that TrackingData and neural time are of the same length
    
    % Input:
    % TrackingData - Structure containing tracking data
    
    % Output:
    % TrackingData - Modified tracking data with equal length
    
    % Check if the length of TrackingData is greater than the length of neural time
    if length(TrackingData.XPosition) > length(TrackingData.neural_time)
        % Get field names of the TrackingData structure
        fields = fieldnames(TrackingData);
        
        % Loop through the fields and adjust their lengths if not a structure
        for i = 1:length(fields)
            if isstruct(TrackingData.(char(fields(i))))
                continue
            else
                % Adjust the length of the field to match neural_time
                TrackingData.(char(fields(i))) = TrackingData.(char(fields(i)))(1:length(TrackingData.neural_time));
            end
        end
    end
end
