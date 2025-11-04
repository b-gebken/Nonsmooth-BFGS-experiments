% The function (8.4) from [Lewis, Wylie (2019)]

function [val,I] = lw2019_84_f(x,g_arr,M_cell,d_arr)

    m = size(g_arr,2);

    val = -Inf;
    tmp1 = g_arr'*x;
    tmp2 = d_arr/24 * norm(x,2)^4;
    for i = 1:m
        tmp_val = tmp1(i) + 1/2 * x'*M_cell{i}*x + tmp2(i);
        if(tmp_val > val)
            val = tmp_val;
            I = i;
        end
    end
    
end

