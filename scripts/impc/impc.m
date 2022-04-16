% This function averages pupil data trials.
% The code is run by using the following command:
% 		impc(<experiment data file>, <sound data file>)
% e.g.
% 		impc('prepared_p1v1b_everything_sample.txt', 'sounds.txt')
%
%     impc - i-modified-pupil-code - just because
% This code will generate three files:
% File 1: Saved processed metadata
% File 2: Uncompiled averaged data
% File 3: Compiled averaged data (data compiling per sentence)
%
% This code will also generate four plots:
% Plot 1: Non-normalized, not blink pruned data per trial 
% Plot 2: Normalized, blink pruned data per trial
% Plot 3: Blinks per trial
% Plot 4: Compiled data per sentence
% 
% This code expects the following attributes in the specified order as input:
% 1. Recording label 
% 2. Right in blink 
% 3. Right in saccade 
% 4. Right pupil diameter
% 5. Timestamp 
% 6. Trial 
% 7. Sound 
% 8. Sentence
%
% Be sure to invoke this script only on previously prepared experiment data files.  
% See the prepare_file.m script for explanations.
%
% You may also want to make use of the prepare_sound_file.m file, which 
% collects sounds from an experiment data file and creates a list of them with 
% predefined start times.
function impc( exp_file_name, sound_file_name )

	% check the number of input arguments, there should be at least two
  if nargin<2
    disp('Usage: impc(<experiment data file>, <sound data file>)');
    error('Missing input! Please provide the experiment data and sound data filenames.');
  end

	% read metadata about the input data file 
	%   pathstr - not used
	%   name - the name of the file itself, it's used to create a directory and 
	% 	      serves as a prefix for evaluation files
	%   ext - not used
	[pathstr, name, ext] = fileparts(exp_file_name);

	% ==============================================
	% CONSTANTS
	% ==============================================

	% time difference between two measurement points in milliseconds
	TIMEDIF = 1;  % in milliseconds
	% the number of points that interest us after a sound is played
	NPOINTS = 1200;
	% relevant time = product of the number of points and the point time resolution
	TIMEREL = TIMEDIF*NPOINTS;

  % disable warnings
  warning('off');
  % disable paging of the output
  more off;
  
  % get the sound data
  sound_data = extract_sound_data(sound_file_name);
  
  % extract experiment data
	%   trial_matrix - the pupil size, a matrix, per trial in a column
	%   blink_matrix - the blinks, a matrix, blink data per trial in a column 
  %   trialnum - number of trials
  %   trial_name_arr - the trial names, by trial index
	%   trial_sentence_arr - the trial sentences, by trial index
	%   trial_sound_arr - the trial sounds, by trial index
	%   trial_counter_arr - the number of relevant records in a trial, by trial index
  printf('Phase 1. Extracting data, please wait ');
  [trialnum, trial_matrix, blink_matrix, trial_name_arr, trial_sentence_arr, trial_sound_arr, trial_counter_arr] = extract_experiment_data(exp_file_name, sound_data, TIMEREL);
	raw_trial_data = trial_matrix;

  printf('Phase 2. Procesing data ... \n');
  
  % split trial data into accepted and rejected trials
  % [accepted_trials,rejected_trials] = split_by_acceptance(trialnum,trial_matrix);
  
  % interpolate data
  trial_matrix = interpolate(trialnum, trial_matrix);
  
	% normalize pupil size variation
	grandmean = mean(mean(trial_matrix)); 
	trial_matrix = trial_matrix/grandmean;
  
  % smooth data
  trial_matrix = smooth(trialnum, trial_matrix);
  
  printf('Phase 3 Starting output ... \n');
  
	% ================================================
	% PLOT THE NORMALIZED DATA
	% ================================================

	figure(1)
	for i=1:trialnum
    	plot(1:trial_counter_arr(i),raw_trial_data(:,i))
    	hold all
	end

	figure(2)
	for i=1:trialnum
    	plot(1:trial_counter_arr(i),trial_matrix(:,i))
    	hold all
	end

	figure(3)
	for i=1:trialnum
    	plot(1:trial_counter_arr(i),blink_matrix(:,i))
    	hold all
	end

	% ================================================
	% SAVE PROCESSED DATA
	% ================================================

	% save trial metadata
	mkdir(name);
	trialfile = [name '\trials.xls']
	fid=fopen(trialfile, 'wt');
	for i=1:trialnum
    	fprintf(fid,'%d\t%s\t%s\t%d\t%s\n',i, trial_name_arr{i},trial_sentence_arr{i},trial_counter_arr(i),trial_sound_arr{i});
	end
	fclose(fid);

	% save uncompiled data
	uncompiledfile = [name '\uncompiled.xls']
	fid=fopen(uncompiledfile, 'wt');
	for i=1:trialnum
    	fprintf(fid,'%s',trial_sentence_arr{i});
    	fprintf(fid,'\t');
	end
	fprintf(fid,'\n');

	for i=1:trialnum
    	fprintf(fid,'%s',trial_sound_arr{i});
    	fprintf(fid,'\t');
	end
	fprintf(fid,'\n');

	maxfilecount = max(trial_counter_arr);     
	for j=1:maxfilecount
    	for i=1:trialnum
        	if j<=trial_counter_arr(i)
            	fprintf(fid,'%.6f',trial_matrix(j,i));
        	else
            	fprintf(fid,'%.6f',0.0);  % print 0 if beyond max number of averaged values
        	end        
        	fprintf(fid,'\t');
    	end
    	fprintf(fid,'\n');
	end
	fclose(fid);

	% save compiled data
	% make list of unique sentences and add similar data together averaged over
	% all occurences
	compiledfile = [name '\compiled.xls']
	scount=1;
	unisentences{scount}=trial_sentence_arr{1};
	sencount(scount) = 1;
	for i=2:trialnum
    	found=0; 
    	j=1;
    	while ((found==0) && (j<=scount)) 
        	if strcmp(trial_sentence_arr{i},unisentences{j})
            	found=1;
            	senindex = j;
        	end
        	j=j+1;
    	end
    	if found==1
        	encount(senindex) = sencount(senindex)+1;
    	else
        	scount = scount+1;
        	unisentences{scount}=trial_sentence_arr{i};
        	sencount(scount)=1;
    	end
	end
	
	%now compile the data
	comparr = zeros(NPOINTS,scount);
	for i=1:scount
    	tcount(i)=NPOINTS; %max points in compiled data
    	for j=1:trialnum
       		if strcmp(trial_sentence_arr{j},unisentences{i}) 
          		if trial_counter_arr(j)<tcount(i)
              		tcount(i) = trial_counter_arr(j); % update max points
          		end
          		for k=1:tcount(i)
              		comparr(k,i)=comparr(k,i)+trial_matrix(k,j);
          		end
       		end
    	end
    	comparr(:,i)=comparr(:,i)/sencount(i);  %average over all occurrences of this sentence
	end

	% plot compiled data
	figure(4)
	for i=1:scount
    	plot(1:tcount(i),comparr(:,i))
    	hold all
	end

	% save the compiled data
	fid=fopen(compiledfile, 'wt');
	for i=1:scount
    	fprintf(fid,'%s',unisentences{i});
		fprintf(fid,'\t');
	end
	fprintf(fid,'\n');
	maxfilecount = max(tcount);     
	for j=1:maxfilecount
 		for i=1:scount
    		if j<=tcount(i) 
        		fprintf(fid,'%.6f',comparr(j,i));
    		else
        		fprintf(fid,'%.6f',0.0);  % print 0 if beyond max number of avergaed values
    		end        
    		fprintf(fid,'\t');
 		end
 		fprintf(fid,'\n');
	end
	fclose(fid);
end




% Extracts the experiment data from the data file, returns it in the form of arrays and matrices.
% @param experiment_file the file containing the experiment data
% @param sound_data the previously extracted sound data with sound begin times
% @param timerel relevant time span to be examined (in milliseconds)
% @return the experiment data in the form or arrays and matrices
%   trialnum - total number of trial in experiment file
%   trial_matrix - pupil size matrix (columns represent trials, rows points of time in trial)
%   blink_matrix - blink matrix (columns represent trials, rows points of time in trial)
%   trial_name_arr - array of trial names indexed by trial index
%   trial_sentence_arr - array of trial sentences indexed by trial index
%   trial_sound_arr - array of trial sounds indexed by trial index
%   trial_counter_arr - the number of relevant records in a trial
function [trialnum,trial_matrix,blink_matrix,trial_name_arr,trial_sentence_arr,trial_sound_arr,trial_counter_arr] = extract_experiment_data( experiment_file, sound_data, timerel )
  % input format 
  I_FORMAT_SPEC = '%s %d %d %s %f %s %s %s';
  % segment size when reading huge data files
	SEGMENT_SIZE = 1000;
  
  % open the experiment data
 	fid = fopen(experiment_file);
  	
  % return variables
  trialnum = 0;
  trial_matrix = [];
  blink_matrix = [];
  trial_name_arr = [];
  trial_sentence_arr = [];
  trial_sound_arr = [];
  trial_counter_arr = [];
  
  % local variables
	soundref = '';
	begin_timerel = 0;
	starttime = 0;
	counter = 0;
  pupilrec = [];
  blinkrec = [];

	while ~feof(fid)
		expdata = textscan(fid, I_FORMAT_SPEC, SEGMENT_SIZE, 'delimiter','\t');
		segment_length = length(expdata{1});

    if segment_length>0
      rec_idx = 1;
      [recording,blink,saccade,pupil,timestamp,sentence,sound,trial] = get_record(rec_idx,expdata);

      while rec_idx<=segment_length
        if (not(strcmp(sound, soundref)))
          if ( counter > 0 )
            % store current pupil data
            trial_counter_arr(trialnum) = counter;
            trial_matrix(:,trialnum) = pupilrec;
            blink_matrix(:,trialnum) = blinkrec;
            % reset data
            counter = 0;
            pupilrec = [];
            blinkrec = [];  
          end 
          
          % find begin of relevant time frame
          n=strmatch(sound,char(sound_data{1}),'exact');
          if (length(n) > 0 ) 
            soundref = char(sound_data{1}(n));  
            begin_timerel = sound_data{2}(n);
          else
            printf('Warning! Missing begin time for sound [%s] in sound data file! Start time will be set to 0 ms.', sound);
            soundref = sound;
            begin_timerel = 0;
          end
          
          starttime = timestamp;    
          trialnum = trialnum+1;
        end

        % trial meta data
        trial_name_arr{trialnum} = trial;
        trial_sentence_arr{trialnum} = sentence;
        trial_sound_arr{trialnum} = sound;
     
        endtime = begin_timerel+timerel;
        while ((rec_idx<=segment_length) && (strcmp(sound,soundref)))
          rectime = timestamp-starttime;
          
          %printf('.... Recordtime: %d => [ Starttime: %d  | Endtime: %d ] .... \n', rectime, begin_timerel, endtime);
          if ((rectime>=begin_timerel) && (rectime<endtime))  
            counter = counter+1;
            pupilrec(counter) = pupil;
            blinkrec(counter) = blink;
          end
          
          rec_idx = rec_idx+1;
          if rec_idx<=segment_length
            [recording,blink,saccade,pupil,timestamp,sentence,sound,trial] = get_record(rec_idx,expdata);
          end
        end
      end

      % give some progress feedback
      printf('.');
    end
	end
  
  if ( counter > 0 )
    % store last pupil data
    trial_counter_arr(trialnum) = counter;
    trial_matrix(:,trialnum) = pupilrec;
    blink_matrix(:,trialnum) = blinkrec;
  end

 	fclose(fid);  
  printf('\n');
end




% Extracts the sound data, sound begin times, from the sound file.
% @param sound_file_name the name of the file containing sound data
% @return an array of attributes
%   sound_data - the sound data and parameters
function sound_data = extract_sound_data( sound_file_name )
	fid = fopen(sound_file_name);
	sound_data = textscan(fid, '%s %d', 'delimiter','\t');
	fclose(fid);
end




% Splits the trials into rejected and accepted trials depending on whether the number or values deviating 
% from the mean by more than 3 standard deviations are more or less than 15%.
% @param trialnum the number of trials
% @param trial_data the trial data; pupil size by trial number in each column
% @return two indexed of accepted and rejectd trials
%   accepted_trials - the indices of the accepted trials
%   rejected_trials - the indices of the rejected trials
function [accepted_trials,rejected_trials] = split_by_acceptance( trialnum, trial_data )
  accepted_trials = [];
	rejected_trials = [];
  
	for trial_id = 1:trialnum
		% find the mean of the current trial
      me = mean(trial_data(:,trial_id));
		% find the standard deviation of the current trial
      stdev = std(trial_data(:,trial_id));
		% find the points in the trial that are less than 3 std's
      rejected_rows = find(trial_data(:,trial_id) < (me - 3*stdev)); 
    
		% set to zero all the trials that deviate from 3 std from the mean
    	trial_data(rejected_rows,trial_id) = 0;
		% check trial - whether we accept of reject
    	if (0.15*length(trial_data(:,trial_id)))>length(rejected_rows) 
        	accepted_trials = [accepted_trials trial_id];
    	else
        	rejected_trials = [rejected_trials trial_id];
    	end
	end
end




% Interpolats the vectors in the data matrix.
% @param trialnum the number of trials
% @param trial_matrix the trial data matrix
% @return the interpolated matrix data
function interpolated_matrix = interpolate( trialnum, trial_matrix )
  for trial_id = 1:trialnum
    trial_data = trial_matrix(:,trial_id);
    trial_matrix(:,trial_id) = doInterpolate(trial_data);
  end
  interpolated_matrix = trial_matrix;
end




% Interpolates a single data vector.
% @param input_data the input data
% @return interpolated_data the interpolated data
function interpolated_data = doInterpolate( input_data )
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
  
  if( length(skips) > 2 )
    even_skips = 0;
  end

  % loop through skips;-1 is necessary because the end of all skips is registered as a skip
	for tau = 1:length(skips)-1-even_skips,
    	% fill the points -66 : 00000 : 132 with the interpolated line ->
    	% using the equation y-y0 = m(x-x1)
      zeros_trailing_range_end = zeros_trailing(skips(tau+1)-1)+samp2+samp1-1;
      zeros_trailing_range_begin = zeros_trailing(skips(tau))-1;
      
      zeros_range_end = zeros(skips(tau+1)-1)+samp2+samp1-1;
      zeros_range_begin = zeros(skips(tau))-1;
      
      i_vec_end_val = extended_data(zeros_trailing_range_end);
      i_vec_beg_val = extended_data(zeros_trailing_range_begin);
      i_vec_ini_val = extended_data(zeros_range_begin);
      
    	interpolated_vector = (((i_vec_end_val-i_vec_end_val) / (zeros_range_end- zeros_range_begin)) * 
        ((zeros_range_begin:zeros_range_end)-zeros_range_begin)) + i_vec_ini_val;
      
      extended_data(zeros_trailing_range_begin:zeros_trailing_range_end) = interpolated_vector';
	end
	interpolated_data = extended_data(samp1:end-samp2);
end




% Smooths the trial data in the trial matrix.
% @param trialnum the trial number
% @param trial_matrix the trial data
% @return the trial matrix with smooth-ed trial data
function smooth_matrix = smooth( trialnum, trial_matrix )
  for trial_id = 1:trialnum
    trial_data = trial_matrix(:,trial_id);
    trial_matrix(:,trial_id) = doSmooth(trial_data);
  end
  smooth_matrix = trial_matrix;
end




% Smooths the interpolated data of the vector. Uses vector convolving to accomplish that.
% @param interpolated_data the interpolated data
% @return the smooth-ed data
function smooth_data = doSmooth( interpolated_data )
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




% Reads a record from a record list and returns an array of values representing the record.
% @param k the record index
% @param expdata all records in the form of a matrix
% @return an array containing the following attributes in the order of their listing:
% 	recording - the recording which was being played
%	  blink - a boolean that indicates whether there was a blink
%	  saccade - movement of the eye
%	  pupil - the pupil size
%	  timestamp - the timepoint of the measurement
%	  sentence - the sentence being played
%	  sound - sound file being played
%	  trial - the trial indicator
function [recording,blink,saccade,pupil,timestamp,sentence,sound,trial] = get_record( k, expdata )
	% the recording is a character string	
	recording = char(expdata{1}(k));
	% the blink is a boolean represented as 1 or 0
	blink = expdata{2}(k);
	% floating point eye direction
	saccade = expdata{3}(k);
	% pupil size is a number, all negative values are eliminated, replaced with a 0
	if length(str2num(char(expdata{4}(k))))>0  % detect .
		pupil = str2num(char(expdata{4}(k)));
	else
		pupil=0;
	end
	% the timestamp is in milliseconds
	timestamp = expdata{5}(k);
	% the sentence is a character string
	sentence = char(expdata{8}(k));
	% the sound is a character string
	sound = char(expdata{7}(k));
	% the trial is a character string
	trial = char(expdata{6}(k));
end