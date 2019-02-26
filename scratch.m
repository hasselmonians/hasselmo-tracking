

%% Repeating Patterns
% Define a vector that repeats n times, such that:
x = repmat([0, 1, 2, 3]', 9, 1);

% then try to model it using the forward difference approach
x0 = [0, 1, 2, 3, 0];

function dx = forward_diff(x, ii, order)
    % evaluates the order'th derivative at index ii of the vector x
    assert(ii > order, 'length of x must be greater than derivative order')

    if order == 0
        dx = x(ii); return;
    elseif order == 1
        dx = -x(ii) + x(ii-1); return;
    elseif order == 2
        dx = x(ii) - 2*x(ii-1) x(ii-2); return;
    elseif order == 3
        dx = -x(ii) + 3*x(ii-1) - 3*x(ii-2) + x(ii-3); return;
    elseif order == 4
        dx = x(ii) - 4*x(ii-1) + 6*x(ii-2) - 4*x(ii-3) + x(ii-4); return;
    end
end

function y = iterate(x, t)
    % given the vector x, predict the rule and iterate forward t steps
    % using the Hille series
    dx1 = forward_diff(x, length(x), 1);
    dx2 = forward_diff(x, length(x), 2);
    dx3 = forward_diff(x, length(x), 3);
    dx4 = forward_diff(x, length(x), 4);

    y = zeros(length(x) + t, 1);
    y(1:length(x)) = x;

end
