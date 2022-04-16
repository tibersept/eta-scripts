% Run as follows:
%   ipolate()
% 
% Runs the interpolation code with some example data
function ipolate()
  % define a sine wave on 1 instead of 0
%  for iterations = 1:4 
    source = [0:1/3000:0.2];
    a = 1;
    f =100;
    points = a*sin(2*pi*f*source)+2;
   
    % gaps in sine wave
    gaps = sort([randi(100),randi(300)+200]);
    
    for gaps_i = 1:2
      gap_beg = gaps(gaps_i);
      gap_end = min(gap_beg + 25,600);
      points(gap_beg:gap_end) = 0;
    end
    
    % plot(source(1:100),points(1:100))
    original_data = points';
    ipolated_data = interpolate(original_data);
    smooth_data = smooth(ipolated_data);
    ipolated_data = ipolated_data';
    smooth_data = smooth_data';
    plot(source(1:600),smooth_data(1:600));
    hold all
    plot(source(1:600),ipolated_data(1:600));
    hold all
%  end
end

% Smoothds the interpolated data
% @param interpolated_data the interpolated data
% @return the smooth-ed data
function smooth_data = smooth( interpolated_data )
	S = 83; % setting the number of smoothing points
	smoo = ones(1,S);
	% performing edge correction
	correlated_data = [flipud(interpolated_data(end-floor(S/2)+1:end)); 
        interpolated_data;  
        flipud(interpolated_data(1:floor(S/2)))];

  % convolve vectors
	convolved_data = conv(correlated_data,smoo/S,'same');
	smooth_data = convolved_data(floor(S/2)+1:end - floor(S/2));
end

% Interpolates the data.
% @param input_data the input data
% @return interpolated_data the interpolated data
function interpolated_data = interpolate( input_data )
	% finding the zero points
	zeros = find(input_data==0);
  
	samp1 = 66; % 66 pts behind the zeros              
	samp2 = 132; % 132 pts after the zeros

	extended_data = [flipud(input_data(1:samp1-1));                  
                    input_data; 
                    flipud(input_data(end - samp2+1:end))];
                    

	% trying to find the length of the zero points 
	zeros_trailing = [zeros; 0];
	zeros_leading = [0; zeros];                    
	skips = find((zeros_trailing-zeros_leading)~=1);

	% check if skips are even
	even_skips = mod(length(skips),2);

  % loop through skips
	for tau = 1:length(skips)-even_skips,
    	% fill the points -66 : 00000 : 132 with the interpolated line ->
    	% using the equation y-y0 = m(x-x1)
    	extended_data(
        zeros_trailing(skips(tau))-samp1+samp1-1:
        zeros_trailing(skips(tau+1)-1)+samp2+samp1-1
      ) = (((
              extended_data(zeros_trailing(skips(tau+1)-1)+samp2+samp1-1)-
              extended_data(zeros_trailing(skips(tau))-samp1+samp1-1)
          ) / (zeros(skips(tau+1)-1)+samp2+samp1-1-(zeros(skips(tau))-samp1+samp1-1))
        ) * ((zeros(skips(tau))-samp1+samp1-1:zeros(skips(tau+1)-1)+samp2+samp1-1)-(zeros(skips(tau))-samp1+samp1-1))
      ) + (extended_data(zeros(skips(tau))-samp1+samp1-1));
	end
	interpolated_data = extended_data(samp1:end-samp2);
end