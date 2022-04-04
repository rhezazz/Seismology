function [misfit] = misfit(t_obs, t_cal)
    ls = length(t_obs);
    for j = 1 : ls
        m(j) = ((t_cal(j) - (t_obs(j))))^2;
    end
    misfit = sqrt((1/ls)*sum(m));
end
