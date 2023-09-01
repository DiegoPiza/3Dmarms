function [fitresult, gof] = nakarushtonDB(xData,yData)
    % Define the Naka-Rushton equation for curve fitting. The Naka-Rushton equation is 
    % typically used in psychophysics to describe the relationship between stimulus 
    % intensity and perception intensity.
    % Equation form: r * (x^n / (x^n + k^n)) + b

    % 'ft' describes the type of curve fit we are using - in this case, the Naka-Rushton equation.
    ft = fittype( '(r*(((x.^n)/((x.^n)+(k^n)))))+b', 'independent', 'x', 'dependent', 'y' );

    % Set up fit options for Nonlinear Least Squares method.
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

    % Turn off display of fit progress.
    opts.Display = 'Off';

    % Set the lower boundary for the parameters [r, n, k].
    opts.Lower = [0 0.2 0.4];

    % Set the starting point for the parameters [r, n, k].
    opts.StartPoint = [0.1 0.7 0.5];

    % Set the upper boundary for the parameters [r, n, k].
    opts.Upper = [1 2.5 0.6];

    % Fit the model to the data using the specified fit type and options.
    % 'fitresult' will contain the fitted curve and 'gof' gives the goodness of fit.
    [fitresult, gof] = fit( xData, yData, ft, opts );
end
