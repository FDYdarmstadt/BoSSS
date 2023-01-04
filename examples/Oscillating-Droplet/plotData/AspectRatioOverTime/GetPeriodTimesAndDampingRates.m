function [ periodTimes, dampingRates ] = GetPeriodTimes( times, values, numP )
%GETPERIODTIMES Summary of this function goes here
%   Detailed explanation goes here

numTs = length(times);

times_max = [];         % times at maxima
values_max = [];
times_min = [];         % times at minima
values_min = [];
times_oneNeg = [];      % times at neg value 1 crossings (+ -> -)
times_onePos = [];      % times at pos value 1 crossings (- -> +)

v_max = values(1);
v_min = -10000;
time_maxmin = 0;
search_max = true;                 % if false the local minima is searched for

for ts=2:numTs-1
    cVal = values(ts);
    nVal = values(ts+1);
    cTime = times(ts);
    if search_max
        if cVal > v_max 
            v_max = cVal;
            time_maxmin = cTime;
        elseif nVal < v_max 
            times_max = [times_max, time_maxmin];
            values_max = [values_max, v_max];
            v_min = cVal;
            search_max = false;
        end       
    else
        if cVal < v_min 
            v_min = cVal;
            time_maxmin = cTime;
        elseif nVal > v_min 
            times_min = [times_min, time_maxmin];
            values_min = [values_min, v_min];
            v_max = cVal;
            search_max = true;
        end   
    end
    pVal = values(ts-1);
    dt = times(ts) - times(ts-1);
    dt1 = times(ts-1) + (1 - pVal)*dt/(cVal - pVal);
    if pVal >= 1 && cVal < 1        
        times_oneNeg = [times_oneNeg, dt1];
    elseif pVal < 1 && cVal >= 1
        times_onePos = [times_onePos, dt1];
    end
end

periodTimes = [];
for i=1:numP
    periodTimes = [periodTimes; times_max(i), times_max(i+1) - times_max(i)];
end
for i=1:numP
    periodTimes = [periodTimes; times_min(i), times_min(i+1) - times_min(i)];
end
for i=1:numP
    periodTimes = [periodTimes; times_oneNeg(i), times_oneNeg(i+1) - times_oneNeg(i)];
end
for i=1:numP
    periodTimes = [periodTimes; times_onePos(i), times_onePos(i+1) - times_onePos(i)];
end

dampingRates = [];
for i=1:numP
    dRate = log((values_max(i+1)-1)/(values_max(i)-1)) / (times_max(i) - times_max(i+1));
    dampingRates = [dampingRates; times_max(i), dRate];
end
for i=1:numP
    dRate = log((values_min(i+1)-1)/(values_min(i)-1)) / (times_min(i) - times_min(i+1));
    dampingRates = [dampingRates; times_min(i), dRate];
end

end

