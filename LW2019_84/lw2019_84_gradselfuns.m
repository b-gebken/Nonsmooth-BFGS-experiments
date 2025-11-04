% The gradients of all selection functions of function (8.4) from [Lewis,
% Wylie (2019)]

function grads = lw2019_84_gradselfuns(x,g_arr,M_cell,d_arr)

m = size(g_arr,2);

grads = g_arr + d_arr/6 .* (norm(x,2)^2*x);

for i = 1:m
    grads(:,i) = grads(:,i) + M_cell{i}*x;
end

end

