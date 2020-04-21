function [lag, R_xy, P_xy, n] = cross_corr(x, y, nlags, bad_flag)

% calculate cross correlations
%
% [lag, R_xy, P_xy, n] = cross_corr(x, y, nlags, bad_flag)
%
% input:
%       x,y the data sets
%       nlags are number of lags to process
%       bad_flag is data value that indicates missing or bad data
%
% output:
%       lag is lag from -nlag to +nlag
%       R_xy are covariances.  
%       P_xy are correlations
%       n is array containing the number of data points used to 
%       calculate each R
%       
%       Subroutine computes correlation between a(t) and b(t+lag). A positive 
%       lag therefore means that a(t) precedes b(t). In other words a(t) leads 
%       b(t).

% put in column format
x=x(:);
y=y(:);

% this many data points
N=length(x);

% initialize output 
R_xy = zeros(2*nlags+1,1);
P_xy = zeros(2*nlags+1,1);
n = zeros(2*nlags+1,1);
% check that y is same length
if (length(y)~=N)
    fprintf('x and y different lengths\n')
    return;
end

% find means
if(isnan(bad_flag))
    bad_flag=1e35;
	id = find(isnan(x));
	x(id) = bad_flag+0*id;
	id = find(isnan(y));
	y(id) = bad_flag+0*id;
end
good_id = find(x~=bad_flag);
if(length(good_id)>0)
    mean_x = mean(x(good_id));
else
    fprintf('no data found\n')
    return;
end
good_id = find(y~=bad_flag);
if(length(good_id)>0) 
    mean_y = mean(y(good_id));
else 
    fprintf('no data found\n') 
    return; 
end

% do the lags
cnt = 0;
for l=-nlags:1:nlags,
    cnt = cnt + 1;

% check for neg./pos lag
    if (l<0)
        k=(-1)*l;
        lag2_id = [1:1:(N-k)]';
        lag1_id = lag2_id+k;

    else
        k=l;
        lag1_id = [1:1:(N-k)]';
        lag2_id = lag1_id+k;
    end

% find good data in x series
    good_id = find( (x(lag1_id)~=bad_flag) );
    Ngoodx = length(good_id);

% continue with this lag if ther are data
    if (Ngoodx>0)
        lag1_id = lag1_id(good_id);
        lag2_id = lag2_id(good_id);

% find good data in y-series where x series was good
        good_id = find( (y(lag2_id)~=bad_flag) );
        Ngood = length(good_id);

% continue only if there are data
        if (Ngood>0)
            n(cnt) = Ngood;
            lag1_id = lag1_id(good_id);
            lag2_id = lag2_id(good_id);

% calculate statistics
% dudley's method
if (1)
    mean_1 = mean(x(lag1_id));
    mean_2 = mean(y(lag2_id));
    z1=x(lag1_id)-mean_1;
    z2=y(lag2_id)-mean_2;
end

% nathaniel's method
if (0)
    z1=x(lag1_id)-mean_x;
    z2=y(lag2_id)-mean_y;
end

% get the normalizing variances
            std_1 = sqrt(z1'*z1/Ngood );
            std_2 = sqrt(z2'*z2/Ngood);

% estimate cov. and corr.
            R_xy(cnt) = z1'*z2/Ngood;
            P_xy(cnt) = R_xy(cnt)/(std_1*std_2);

	end % check for good x -data
    end % check for good y -data
end % all lags

% reproduce lags
lag = [-nlags:1:nlags]';