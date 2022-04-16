function ipupil(expfilename,soundfilename,responsefilename)
	dofilter=0;

	if nargin<2
    	error('Please enter the experiment filename and soundfilename');
	end

	[pathstr, name, ext] = fileparts(expfilename);

	timedif = 1;  % in msec
	npoints = 1200; % number of points after soundtime
	grammarlevel = 3;
	conflevel=3;
	timerel = timedif*npoints;

	% read exp file into cellarray
 	fid = fopen(expfilename);
 	expdata = textscan(fid, '%s %d %d %s %d %s %s %s', 'delimiter','\t','HeaderLines',1); 
 	fclose(fid);

 	explen=length(expdata{1})
 

	% read sound file into cellarray
	fid = fopen(soundfilename);
	sounddata = textscan(fid, '%s %d', 'delimiter','\t');
	fclose(fid);
	numsounds = length(sounddata{1})


	rec = 1;
	trialnum = 0;
	trialarr = [];

	%read first record of new sequence
	[recording,blink,saccade,pupil,timestamp,sentence,sound,trial] = getrecord(rec,expdata);

	while rec<explen
  		%find soundtime
    	n=strmatch(sound,char(sounddata{1}),'exact');
    	soundref = char(sounddata{1}(n));
    	soundtime = sounddata{2}(n);
   		starttime = timestamp;
   		counter = 0;
   		trialnum = trialnum+1;
   		%trial
    
   		trialnamearr{trialnum} = trial;
		trialsentencearr{trialnum} = sentence;
		trialsoundarr{trialnum} = sound;
   
		while ((rec<explen) && (strcmp(sound,soundref)))
      		rectime = timestamp-starttime;
      		if ((rectime>=soundtime) && (rectime<soundtime+timerel))  
          		counter = counter+1;
          		pupilrec(counter) = pupil;
          		blinkrec(counter) = blink;
      		end
      		rec = rec+1;
      		[recording,blink,saccade,pupil,timestamp,sentence,sound,trial] = getrecord(rec,expdata);
   		end
   
   		trialcounterarr(trialnum) = counter;
   		trialarr(:,trialnum) = pupilrec;
   		blinkarr(:,trialnum) = blinkrec;
	end

	trialnum

	rawtrialdata = trialarr;

	% remove blinks 
	blinkoffset=50;
	for tr=1:trialnum
 		m=1;   
 		trc = trialcounterarr(tr);
 		while m<trc
   			if blinkarr(m,tr)==1
     			blinkstart = m;
     			blinkcounter = 0;
     			% find end of blink
     			bn = blinkstart;
     			while ((bn<trc) && (blinkarr(bn,tr)==1))
        			blinkcounter = blinkcounter+1;
        			bn = bn+1;
     			end
     
				%     blinkstart
				%     bn
    			if blinkstart<=blinkoffset
        			bos = 1;
    			else
        			bos = blinkstart-blinkoffset;
    			end

    			if bn+blinkoffset>trc
        			boe = trc;
    			else
        			boe = bn+blinkoffset;
    			end

    			if boe==trc  % then we have run over boundary so use preblink value
         			for i=bos:boe
            			trialarr(i,tr) =  trialarr(bos,tr);
         			end
     			elseif bos==1  % use end of blink as value
         			for i=bos:boe
            			trialarr(i,tr) =  trialarr(boe,tr);
         			end
    			else % use linear interpolation between start en end of blink
         			bc = trialarr(bos,tr);
         			br = (trialarr(boe,tr)-trialarr(bos,tr))/(boe-bos);
         			for i=bos:boe
            			trialarr(i,tr)= (i-bos)*br+bc;
         			end
     			end

     			m=bn;

   			end
   			m=m+1;
 		end
	end

	%get mean and divide all values by the mean to normalize pupilsize variation as
	%result of camera distance

	grandmean = mean(mean(trialarr)); 
	trialarr = trialarr/grandmean;

	figure(1)
	for i=1:trialnum
    	plot(1:trialcounterarr(i),rawtrialdata(:,i))
    	hold all
	end

	figure(2)
	for i=1:trialnum
    	plot(1:trialcounterarr(i),trialarr(:,i))
    	hold all
	end

	figure(3)
	for i=1:trialnum
    	plot(1:trialcounterarr(i),blinkarr(:,i))
    	hold all
	end

	mkdir(name);
	trialfile = [name '\trials.xls']
	% file to write result to
	fid=fopen(trialfile, 'wt');
	for i=1:trialnum
    	fprintf(fid,'%d\t%s\t%s\t%d\t%s\n',i, trialnamearr{i},trialsentencearr{i},trialcounterarr(i),trialsoundarr{i});
	end
	fclose(fid);

	uncompiledfile = [name '\uncompiled.xls']
	fid=fopen(uncompiledfile, 'wt');
	for i=1:trialnum
    	fprintf(fid,'%s',trialsentencearr{i});
    	fprintf(fid,'\t');
	end
	fprintf(fid,'\n');

	for i=1:trialnum
    	fprintf(fid,'%s',trialsoundarr{i});
    	fprintf(fid,'\t');
	end
	fprintf(fid,'\n');

	maxfilecount = max(trialcounterarr);     
	for j=1:maxfilecount
    	for i=1:trialnum
        	if j<=trialcounterarr(i)
            	fprintf(fid,'%.6f',trialarr(j,i));
        	else
            	fprintf(fid,'%.6f',0.0);  % print 0 if beyond max number of averaged values
        	end        
        	fprintf(fid,'\t');
    	end
    	fprintf(fid,'\n');
	end
	fclose(fid);

	compiledfile = [name '\compiled.xls']
	%make list of unique sentences and add similar data together averaged over
	%all occurences
	scount=1;
	unisentences{scount}=trialsentencearr{1};
	sencount(scount) = 1;
	for i=2:trialnum
    	found=0; 
    	j=1;
    	while ((found==0) && (j<=scount)) 
        	if strcmp(trialsentencearr{i},unisentences{j})
            	found=1;
            	senindex = j;
        	end
        	j=j+1;
    	end
    	if found==1
        	encount(senindex) = sencount(senindex)+1;
    	else
        	scount = scount+1;
        	unisentences{scount}=trialsentencearr{i};
        	sencount(scount)=1;
    	end
	end

	unisentences
	sencount
	%now compile the data
	comparr = zeros(npoints,scount);
	for i=1:scount
    	tcount(i)=npoints; %max points in compiled data
    	for j=1:trialnum
       		if strcmp(trialsentencearr{j},unisentences{i}) 
          		if trialcounterarr(j)<tcount(i)
              		tcount(i) = trialcounterarr(j); % update max points
          		end
          		for k=1:tcount(i)
              		comparr(k,i)=comparr(k,i)+trialarr(k,j);
          		end
       		end
    	end
    	comparr(:,i)=comparr(:,i)/sencount(i);  %average over all occurences of this sentence
	end

	figure(4)
	for i=1:scount
    	plot(1:tcount(i),comparr(:,i))
    	hold all
	end

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


%
% Reads a record from a record list.
%
function [recording,blink,saccade,pupil,timestamp,sentence,sound,trial] = getrecord(k,expdata)
         recording = char(expdata{1}(k));
         blink = expdata{2}(k);
         saccade = expdata{3}(k);
         if length(str2num(char(expdata{4}(k))))>0  % detect .
             pupil = str2num(char(expdata{4}(k)));
         else
             pupil=0;
         end
         timestamp = expdata{5}(k);
         sentence = char(expdata{6}(k));
         sound = char(expdata{7}(k));
         trial = char(expdata{8}(k));
end
