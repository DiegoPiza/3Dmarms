function [trainedClassifier, validationAccuracy] = trainClassifierfloor(trainingData, responseData)
    % trainClassifier: Trains a classifier and returns a trained classifier and its accuracy.
    % Input:
    %   trainingData: Matrix of predictor data.
    %   responseData: Vector of response data.
    % Output:
    %   trainedClassifier: A struct containing the trained classifier.
    %   validationAccuracy: Accuracy of the trained classifier.

    % Extract predictors and response
    % Convert input to table
    inputTable = array2table(trainingData);

    % Set up predictor names
    predictorNames = inputTable.Properties.VariableNames;

    % Set up predictors and response
    predictors = inputTable;
    response = responseData;

    % Train a classifier
    template = templateSVM(...
        'KernelFunction', 'linear', ...
        'PolynomialOrder', [], ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    classificationSVM = fitcecoc(...
        predictors, ...
        response, ...
        'Learners', template, ...
        'Coding', 'onevsone');

    % Create the result struct with predict function
    trainedClassifier = struct();
    trainedClassifier.ClassificationSVM = classificationSVM;

    % Calculate validation accuracy
    cvmodel = crossval(classificationSVM, 'KFold', 5);
    validationAccuracy = 1 - kfoldLoss(cvmodel, 'LossFun', 'ClassifError');

    % Display validation accuracy
    disp(['Validation Accuracy: ', num2str(validationAccuracy * 100), '%']);
end
