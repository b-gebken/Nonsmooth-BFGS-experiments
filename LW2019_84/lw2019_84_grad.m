% The gradient of function (8.4) from [Lewis, Wylie (2019)]

function grad = lw2019_84_grad(x,g_arr,M_cell,d_arr)
    
    [~,I] = lw2019_84_f(x,g_arr,M_cell,d_arr);

    grad = g_arr(:,I) + M_cell{I}*x + d_arr(I)/24 * 4*norm(x,2)^2 * x;

end

