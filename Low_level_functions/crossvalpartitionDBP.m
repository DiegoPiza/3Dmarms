function [test_indices, train_indices] = crossvalpartitionDBP(N, n_folds)
    % Calculate the size of each fold
    fold_size = floor(N / n_folds);

    % Initialize matrices to store fold indices
    fold_indices = zeros(fold_size, n_folds);

    % Initialize matrices for test and train indices
    test_indices = false(N, n_folds);
    train_indices = true(N, n_folds);

    % Create continuous partitions for cross-validation
    for i = 1:n_folds
        % Calculate indices for the current fold
        fold_start = (i - 1) * fold_size + 1;
        fold_end = i * fold_size;

        % Set test indices for the current fold to true
        test_indices(fold_start:fold_end, i) = true;

        % Set train indices for the current fold to false
        train_indices(fold_start:fold_end, i) = false;
    end
end
