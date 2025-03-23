%% Task 1

1/0-1/0
(1+eps+eps^2 == 1+eps^2+eps)
realmin/2^55
1-eps/2+eps/2
realmax('single')+realmax('double')
double(realmax('single'))+realmax('double')

a=(10+eps+eps^2==10+eps^2+eps);
b=(1/0)*(2/Inf);
c=(eps*eps^2==eps^3); d = (NaN==NaN)

%% Task 3

function result = myfloats(ekth, mantis, fp_type)
    % Define exponent and mantissa sizes for different floating-point types
    fp_formats = struct('fp64', [11, 52], 'fp32', [8, 23], ...
                        'fp16', [5, 10], 'bfloat16', [8, 7]);

    % Check if the provided floating-point type is valid
    if ~isfield(fp_formats, fp_type)
        fprintf('Error: Invalid floating-point type.\n');
        return;
    end

    % Get the exponent and mantissa sizes for the given type
    exp_size = fp_formats.(fp_type)(1);
    man_size = fp_formats.(fp_type)(2);

    % Validate input exponent and mantissa length
    if ekth < -(2^(exp_size-1) - 1) || ekth > (2^(exp_size-1) - 1)
        fprintf('Error: Exponent out of range for %s.\n', fp_type);
        return;
    end
    if length(mantis) ~= man_size
        fprintf('Error: Mantissa length does not match %s format.\n', fp_type);
        return;
    end

    % Convert mantissa binary vector to decimal fraction
    mantissa_decimal = sum(mantis .* 2.^(-1:-1:-man_size));

    % Compute the final floating-point value
    result = (1 + mantissa_decimal) * 2^ekth;
    
    % Display result
    fprintf('Floating-point value in decimal: %e\n', result);
end

ekth = 3;  % Exponent value
mantis = [1 0 1 0 0 1 0 1 1 0 0 1 1 0 1 0 1 0 0 1 0 0 1];  % 23-bit mantissa
fp_type = 'fp32';  % Floating-point type

% Call the function
myfloats(ekth, mantis, fp_type);

%% Task 5

syms x y z
f = x^2 + x*y + 2*z;
J = jacobian(f, [x, y, z]);
disp(J);

%% Task 6

syms x y z
f = [2*x + 3*y - 1; 3*x - 2*y - 2];
J = jacobian(f, [x, y, z]);
disp(J);

%% Task 7

format long
x = 1e20;
y = 10;
z = -1e20;

% Standard computation
result = (x + y) + z;
disp(result);

% Small perturbation to see the exact computation
x_pert = 1e20 + eps(1e20); % Slightly different value
y_pert = 10;
z_pert = -1e20;

result_pert = (x_pert + y_pert) + z_pert;
disp(result_pert);

%% Task 9

function P_approx = horner_eval(a, x)
    % Horner's method for polynomial evaluation
    % a: Coefficients of the polynomial [a_n, ..., a_1, a_0]
    % x: Point at which to evaluate the polynomial
    
    n = length(a);
    P_approx = a(1); % Start with the highest-degree coefficient
    
    for i = 2:n
        P_approx = P_approx * x + a(i);
    end
end

% Example usage
a = [1, -3, 2, 4]; % Polynomial: x^3 - 3x^2 + 2x + 4
x = 1.5; % Evaluate at x = 1.5

P_exact = polyval(a, x); % Using MATLAB's built-in function for comparison
P_horner = horner_eval(a, x);

disp(['Horner’s method result: ', num2str(P_horner)]);
disp(['MATLAB polyval result: ', num2str(P_exact)]);
disp(['Error: ', num2str(abs(P_horner - P_exact))]);

%% Task 10

format long;

% Define original values
a = 1.0; b = 2.0; c = 3.0;
d = 1e-8; e = 2.0; f = 3.0; % Small d can cause instability

% Compute the original matrix
F_original = [a*d, a*e, a*f;
              b*d, b*e, b*f;
              c*d, c*e, c*f];

% Introduce a small perturbation in d
perturbation = 1e-10; 
d_perturbed = d + perturbation;

% Compute the new matrix with perturbed d
F_perturbed = [a*d_perturbed, a*e, a*f;
               b*d_perturbed, b*e, b*f;
               c*d_perturbed, c*e, c*f];

% Compute relative error
relative_error = norm(F_perturbed - F_original) / norm(F_original);

% Display results
disp('Original matrix:');
disp(F_original);

disp('Perturbed matrix:');
disp(F_perturbed);

disp(['Relative error in output: ', num2str(relative_error)]);

%% Task 11

format long;

% Define values of x = 10^p for p = 1 to 40
p_values = 1:40;
x_values = 10.^(-p_values);

% Compute e^x - 1 using different methods
direct_exp = exp(x_values) - 1;   % Direct computation
stable_expm1 = expm1(x_values);   % MATLAB’s stable function
taylor_approx = x_values + (x_values.^2)/2;  % First two terms of Taylor series

% Display results for selected values
disp('   p       exp(x)-1          expm1(x)        Taylor Approximation');
for i = 1:length(p_values)
    fprintf('%4d  %14.12e  %14.12e  %14.12e\n', p_values(i), direct_exp(i), stable_expm1(i), taylor_approx(i));
end