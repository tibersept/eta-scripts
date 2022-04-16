% This function averages pupil data trials according to user defined conditions.
%
% The code is run by :
% output = analyse_pupil('input_filename',trial_length,trial_ID)
%
% where input_filename has to be a '.mat' file
% the trial_length is in samples (and NOT in seconds)
% the trial_ID is defined by including the TTL_ID between square brackets, 
% each separated by space  
%
% e.g
% p = analyse_pupil('temp.mat',6001,[101 102 103 104 105]);
%
% This code will generate three files:
% 1. An .xls file which contain the average of the accepted interpolated trials
% 2. An .xls file which contain the average of the accepted non-interpolated trials
% 3. An .txt file which contain the number of accepted and rejected trials

% INFO ::: data input structure MUST be in the following format:
% data_structure consisted of 7 cols:
% 1: Subject No
% 2: Time Stamp
% 3: Pupil Diameter
% 4: Trial No
% 5: TTL
% 6: Right Blink/Saccade
% 7: Right Blink/Saccade
function a = analyse_pupil(filename,trial_size,ID)

	% loads a matlab file => will contain the input data as a matrix
	s = load(filename);
	% creates a string s.xyz 
	oooo = cat(2,'s.',filename(1:end-4));

	% initialisation of constants
	test_orig = eval(oooo);
	trialNO_col = 4;
	data_col = 3;
	test = test_orig;
	accept = [];
	reject = [];

	% for each trial do the following for the no. of trials
	% Split the trials into rejected and accepted trials
	% ... depending on whether the pupil diameter in any case
	% ... has a deviation of more than 3 standard deviations
	for trials = 1:test(end,trialNO_col),
		% get the trial that will be analysed, i is a set of indexes 
    	i = find(test(:,trialNO_col)==trials); 
		% find the mean of the current trial
    	me = mean(test(i,data_col)); 
		% find the standard deviation of the current trial
    	stdev = std(test(i,data_col)); 
		% find the points in the trial that is less than 3std
    	i1 = find(test(i,data_col)< (me - 3*stdev)); 
    	% i2 = find(test(i,5)> (me + 3*stdev));
    	i_tot = [i1+i(1)-1];% i2+i(1)-1];
    
    	% debugging code : trying to identify the length of the zero values
    	% length(i_tot);
    
		% set to zero all the trials that deviates from 3 std from the mean
    	test(i_tot,data_col) = 0; 
		% check trial - whether we accept of reject
    	if 0.15*length(i)>length(i_tot), 
        	accept = [accept trials];
    	else
        	reject = [reject trials];
    	end
	end


	% figure(1) % plotting the cont data before and after the nullifying process	
	% plot(test_orig(:,data_col),'r');hold on
	% plot(test(:,data_col),'g');
	% hold off;

	%%%%%%%%%%%%%%%%%%%%%
	%Interpolation Code %
	%%%%%%%%%%%%%%%%%%%%%

	% get the nullified data
	t = test(:,data_col);
	t1 = t;
	% finding the zero points
	i = find(t==0);
	samp1 = 66; % reference to the 66 pts behind the zeros
	samp2 = 132; % reference to the 132 pts after the zeros


	t_2 = [flipud(t(1:samp1-1)); t; flipud(t(end - samp2+1:end))];

	% trying to find the length of the zero points 
	a1 = [i; 0];
	a2 = [0; i];
	ind = find((a1-a2)~=1);

	% trying to see whether the zeros start and end points are even as
	% otherwise process only up the last matching start and end points
	a = mod(length(ind),2);


	for tau = 1:length(ind)-1-a, %loop for the number of matched start and end points
    	%filling the points -66 : 00000 : 132 with the interpolated line ->
    	%using equation y-y0 = m(x-x1)
    	t_2(a1(ind(tau))-samp1+samp1-1:a1(ind(tau+1)-1)+samp2+samp1-1) = (t_2(a1(ind(tau+1)-1)+samp2+samp1-1)-t_2(a1(ind(tau))-samp1+samp1-1))/(i(ind(tau+1)-1)+samp2+samp1-1-(i(ind(tau))-samp1+samp1-1))*((i(ind(tau))-samp1+samp1-1:i(ind(tau+1)-1)+samp2+samp1-1)-(i(ind(tau))-samp1+samp1-1))+t_2(i(ind(tau))-samp1+samp1-1);
	end
	t1 = t_2(samp1:end-samp2);



	% %Uncomments below if you need to plot before and after the linear
	% %interpolations

	% figure(2) % Compare the between the interpolated cont data and the nullified data
	% plot(t);hold on;
	% plot(t1,'r')
	% hold off;

	%%%%%%%%%%%%%%
	% Smoothing
	%%%%%%%%%%%%%%

	S = 83; % setting the number of smoothing points
	smoo = ones(1,S);
	% performing edge correction
	t2 = [flipud(t1(end-floor(S/2)+1:end)); t1;  flipud(t1(1:floor(S/2)))];

	t3 = conv(t2,smoo/S,'same');
	res = t3(floor(S/2)+1:end - floor(S/2));

	% figure(4) %Plotting between the smoothed interp and non-smoothed interp cont. data
	% plot(t1,'r');hold on
	% plot(res);

	%%%%%%%%%%%%%%%%%
	% Averaging
	%%%%%%%%%%%%%%%%%%

	%Setting up variables to store the trial according to their corresponding
	%TTL ID

	% For interp
	for u = 1:length(ID),
    	pupil{u} = []; %non_interp
    	pupil_np{u} = []; %interp
	end

	%Getting the non-interp data
	res0 = test(:,data_col);

	ID = [101 102]; % those are the only conditions of interest

	for op = 1:length(accept),%do for the number of accepted trials
    	i_acc = find(test(:,trialNO_col)==accept(op)); % get the index of the corresponding accepted trials
    	c = test(i_acc(1),5); % get the ttl code
    	if length(i_acc) == trial_size, % check if the length of the trial has a size of the defined length
        	for uin = 1:length(ID),
            	if c == ID(uin)
                	pupil{uin} = [pupil{uin} res0(i_acc)];
                	pupil_np{uin} = [pupil_np{uin} res(i_acc)];
            	end
        	end
    	else
        	% proceed to the next accepted trial if the accepted trial length
        	% is not the same as the trial_size defined by the user
    	end
	end

	%calculating the mean of the interp (calculates the row means and stores them in me_ave_pupil)
	for yo = 1:length(ID)
    	me_ave_pupil{yo} = mean(pupil{yo}',1);
    	me_ave_pupil0{yo} = mean(pupil_np{yo}',1);
	end

  % index vecotr for number of trials for first condition (TTL) 
  tim = 1:length(me_ave_pupil{1});
  
  %apparently the test_orig matrix contains condition ids or something similar at position 2,1
  %the following code generates names like sub<whatever is at position 2,1>_pupil_interp 
	nam = cat(2,'sub',num2str(test_orig(2,1)),'_pupil_noninterp'); 
	nam0 = cat(2,'sub',num2str(test_orig(2,1)),'_pupil_interp');

	temp_fin = tim';
	temp_fin0 = temp_fin;

	%averaging for the no. of conditions; 
  %creates an indexed matrix of row means; the matrix -> [idx means_101 means_102 ...]
  %assumes row means have the same length, i.e. there are equal number of trials for all conditions 
	for taut = 1:length(ID)
    	temp_fin = [temp_fin me_ave_pupil{taut}']; 
    	temp_fin0 = [temp_fin0 me_ave_pupil0{taut}'];
	end

  %adds additional rows to the bottom of the output matrix; 
  %first additional row contains just 0 to mark the sepration of real data from the max values
  %... which max values (of the rows, i.e the max mean) are added in the second row, the index for that row is 0
	temp_fin = [temp_fin;zeros(1,length(ID)+1); 0 max(temp_fin(:,2:end),[],1)];
	temp_fin0 = [temp_fin0;zeros(1,length(ID)+1); 0 max(temp_fin0(:,2:end),[],1)];

	%writing the results into ".xls"
	xlswrite(nam, temp_fin);
	xlswrite(nam0, temp_fin0);

	nam_acc_rej = cat(2,'sub',num2str(test_orig(2,1)),'_accept_reject_info.txt');

	fileID = fopen(nam_acc_rej,'w+');
	fprintf(fileID,'accepted trials : %d\nrejected trials : %d',length(accept),length(reject));
	fclose(fileID);

	%% CODE FOR SAVING INTO CONT FORMAT
	% out = test_orig;
	% out(:,data_col) = res;
	% 
	% name = cat(2,'cont_data_sub',num2str(durationAF(2,1)),'.txt');
	% 
	% save(name, 'out','-ASCII'); 

	a='OK';
