Sample data in .mat format, loads a workspace in MATLAB with the following variables:
- TrackingData= Data structure with raw tracking data ( XYZ Position (in
meters) and head orientation in rotation quaternions), the neuraltime field is in
seconds (pulse used to sync tracking to neural data, at 60 Hz frequency, each pulse indicates the corresponding blackrock cerebus time for each tracking frame )
- gap= sampled at 30khz, logical array (true = gap in neural data) gaps are
dropped signal from the headstage.
- units= has single and multi units, timestamps are in seconds (blackrock cerebus time).

