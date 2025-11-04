% A wrapper for the structs problem_data and algo_options to convert the
% inputs to symbolic vpa values.

function [problem_data,algo_options] = vpa_wrapper(problem_data,algo_options,vpa_digits)

    digits(vpa_digits);
    problem_data.x0 = vpa(problem_data.x0);
    algo_options.H0 = vpa(algo_options.H0);
    algo_options.descent_threshold = vpa(algo_options.descent_threshold);
    algo_options.step_threshold = vpa(algo_options.step_threshold);
    algo_options.c1 = vpa(algo_options.c1);
    algo_options.c2 = vpa(algo_options.c2);

end

