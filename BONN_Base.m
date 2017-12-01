
%% 
% make sure that you have this algorithm and any2csv.m in the working file 
% folder containing the data files
%%

clear all
clear global               % clear all variables

tic

files = dir('*.mat');  
disp('files in pathway:       ')  
disp(numel(files))
%    
%for a = 1:numel(files) 
%        readout = {a files(a).name};
%        disp(readout)
%end
%
%g = input ('files to analyze  '); 
    % specify files based on readout shown. 
    % to select sequential files, format x:y 
    % -or- put in a array like [a b c e]
g = 1:500;
h = 301:500;

%h = input ('baseline file:    ');
    % select baseline file(s)
    % if multiple baselines, format x:y (sequential) or [a b d]
    
%for manual input
%f = input ('sampling frequency:    ');
    %    f = 400;       % KA sampling frequency for data = 400 Hz
%CH = input('Channels:   ');
    %    CH = 1;        % KA channels = 1
%win0 = input ('sampling window:    ');       
    %    win0 = 100;    % KA window size = 100, eqivalent to 250 ms
%DC = input ('Decomp Level:   ');             
    %    DC = 4;        % KA decomp windo = 4, equivalent to 12.5-25 Hz

    %For UW Madison Data: 
    f = 173.61; %sampling rate is 1024 Hz
    CH = 1; %update as needed
    win0 = ceil(f/4); %250 ms = 1024/4
    DC = 4; %need to also try 4?
    
    %for Mayo Mouse Data:
   % f = 400; %sampling rate is 1024 Hz
   % CH = 1; %update as needed
   % win0 = 100; %250 ms = 1024/4
   % DC = 4; %decomposition to level 4, approx 12-15 hz?
   
    
    % adjust window size so it can be evenly divided by the downsampling
    % factor identified by the decomposition level
    dc = 2^DC;
    win = ceil(win0/dc)*dc;
    w = win/dc;
    
%ww=1; % is used to set up the Decomp Variable - WTF do I need to set this for?


% Set up variables for finding the threshold, h is one number, so length=1.
% BaseStatMedian = zeros(1,length(h));
% BaseStatStDev = zeros(1,length(h));

% clear win0


 
% Find median, standard deviation of baseline files. These values will be 
% used to set the threshold for event identification. 
% 
% 
% 
% Perform this loop on baseline file(s) only 
% 
% h = baseline file id (from input above)

%%
for k = h   
    
%%
        load (files(k).name); % load the k/h-th (baseline) file in files 
       % EEG = clean; % clean is the EEG variable, all channels.
         j = k-h(1)+1;
  %%
 
  
        for col = 1:CH    
                    % perform this loop on each channel individually
   %%     
            [C,L]  = wavedec(EEG(:,col),DC,'db4');
                    % perform wavelet decomposition on the designated channel
                    % to the indicated decomposition level (DC) using the 
                    % db4 wavelet
                    % for DC = 4 on 400 hz data, pull out XXXXX frequency
                    % for DC = 5 on 1024 hz data, pull out yyyyy frequency
                    
                    % From help wavedec
                    % The output vector, C, contains the wavelet decomposition. 
                    % L contains the number of coefficients by level.
                    % C and L are organized as:
                    % C      = [app. coef.(N)|det. coef.(N)|... |det. coef.(1)]
                    % L(1)   = length of app. coef.(N)
                    % L(i)   = length of det. coef.(N-i+2) for i = 2,...,N+1
                    % L(N+2) = length(X).
               
               
            %L2 = [0; L]; % column of L with 0 added to the beginning
         %   Decomp = C((sum(L2(1:ww-1))+1):L(ww)); 
         Decomp = C(1:L(1));
                [rows,~] = size (Decomp); % '~' voids the 1 column. 
                    % Set Decomp equal to the approximation at the DC level
                    % The first value of L gives the number of coefficients 
                    % from the beginning of C that contains the 4th level
                    % approximation from the transform. This contains the 
                    % "most information" while reducing noise.

%
                   
                    
                    % find distance between datapoints n and n-1
        Dist = zeros(rows,1); 
            % creates a column vectors full of zeros that is the size of 
            % the decomposed/approximation of the EEG (rows).
            for m = 2:(rows-w) % size of Decomp (rows) - window (win = 96)
                %w = win/16 or win/DC^4
                
                Dist((m-1)) = abs(Decomp(m)-Decomp(m-1));
                    % modifies Dist and inputs the individual line lengths. 
                    % Starts from the 1st row in Dist and gives the
                    % absolute distance between the 2nd point and 1st point
                    % of the approximation/decomposed EEG.
                    % Ends at rows-6-1 and gives the abs diff between
                    % rows-6 and rows-6-1
            end
     %                              
                    % sum distance between datapoints across window
            LLength = zeros(length(Dist),1); 
                % creates LLength full of a column of zeros that has a
                % length of Dist which is the size of the decomp EEG
            LLStart = sum(Dist(1:w,1)); 
                % sums up line lengths of first window (6 numbers long,
                % represents ~100 datapoints from window)
                % which is Dist 1 to Dist 6,
                % which are line lengths between points 1 to 97 of the 
                % decomp EEG
            LLength(1,1) = LLStart;
                % modifies the first entry of LLength to the sum of line
                % lengths of first window
  %              
        for y = 2:rows-w 
            % this is the new LLength entry row, aka the 2nd window
                % which is sum of Dist 2 to Dist 97
                % which are line lengths between points 2 to 98 of the
                % decomp EEG
            % same as m, the number of times we calculated distance
                % runs from from 2 to rows-window(96) 
                % PROBLEM: the last window, from rows-96 to rows has 0 dist
                % because it wasn't calculated in Dist above.
            yy = y+w-1; 
                % starts at 2+96-1=97 to rows-96+96-1=rows-1
                    % specifies the new Dist entry to add to this 
                    % new window.
            LL_a = LLength(y-1)-Dist(y-1); 
                % LLength(2-1) - Dist(2-1),
                % DLLength(rows-96-1) - Dist(rows-96-1)   
                    % takes the LL sum of previous window and subtracts
                    % the first Dist entry of that previous window
            LLength(y) = LL_a + Dist(yy);
                % Adds the next Dist entry to the new window 
                % Modifies LLength with this new sum of dist in new window
        end
        %
        
    %    j = k-h(1)+1;
            % j = baseline file number in 'files' 
            % - first baseline file number in 'files' 
            % + 1, answer is 1 if there is only one baseline file 
            % used if there are multiple baseline files.
                    
                    % find median and standard deviation of calculated
                    % line lengths
       
            %BaseStatMedian(j,col) = median(LLength(:,col));
            %BaseStatStDev(j,col) = std(LLength(:,col));
            
            BaseStatMedian(j,col) = median(LLength);
            BaseStatStDev(j,col) = std(LLength);
            
            
        end
%%
    clear EEG clean Decomp rows Dist LLength LLStart LL_a y yy C L
 
%%
    
end

%if length(h)>1 
        % the stats are averaged across baseline files and logged into a
 
 % BaseStats variable.
 %%
for b = 1:CH
    BaseStats(1,b) = mean(BaseStatMedian(:,b));
    BaseStats(2,b) = mean(BaseStatStDev(:,b));
end
%else
 %   BaseStats(1,:) = BaseStatMedian;
  %  BaseStats(2,:) = BaseStatStDev;
%end

%clear BaseStatMedian BaseStatStDev j


%%

SF = 3.5;       % SF can be adjusted to increase or decrease the threshold 
               % used to determine if a line length is different from the 
                % median baseline line length. SF is a multiplier and
                % identifies how many standard deviations from the median a
                % line length value must be in order for it to be pulled
                % out in the analysis.
                
%                SF = 1.5;
%%
for k = g
    
    
  %%
    load (files(k).name);
        input_name = files(k).name;            
        [~, name, ~] = fileparts (input_name); 
    %    EEG = clean;  
            fr=length(EEG)/f;
            t=(1/f:1/f:fr);              
            t=t(:); 
    %%
    for col = 1:CH  % run all channels
                    
        [C,L]  = wavedec(EEG(:,col),DC,'db4');
          % perform wavelet decomposition to the indicated 
                    % decomposition level (DC) using the db4 wavelet
            L2 = [0; L];
             Decomp = C(1:L(1));
        %Decomp = C((sum(L2(1:ww-1))+1):L(ww));
            [rows, ~] = size (Decomp);
            % Set Decomp equal to the approximation at the DC level
         
            % find distance between datapoints n and n-1
        Dist = zeros(rows,1);
            for r = 2:rows-w
                Dist(r-1) = abs(Decomp(r)-Decomp(r-1));
            end
            %%
            % sum distance between datapoints across window
            LLength = zeros(rows,1);
            LLStart = sum(Dist(1:w));
            LLength(1) = LLStart;
                for r = 2:rows-w
                    LL_a = LLength(r-1)-Dist(r-1);
                    LLength(r) = LL_a+Dist(r+w);
                end
            LLAve = LLength-BaseStats(1,col);
        
            clear rows LLAveSD
            [rows, ~] = size(Dist);
      
            LLAveSD = zeros(rows,1);
            LLAveSD = LLAve/(SF*BaseStats(2,col)); 
         
        % create binary vector of thresholded data
            [rows, columns] = size(LLAveSD);
                LLTh = zeros(rows,1);
                for r = 1:rows
                    if LLAveSD(r) >= 1
                        LLTh(r) = 1;
                    else, LLTh(r) = 0;
                    end
                end
                
        % erode (open) data to remove noise - half-window size = 1.
            LLTh2 = zeros(size(LLTh));
            nhood = 1;
            GG = length(LLTh);
                for r = 1
                    if sum(LLTh(1:(r+nhood))) == length(LLTh(1:(r+nhood)))
                        LLTh2(r) = 1;    
                    else, LLTh2(r) = 0;
                    end
                end
                for r = 2:GG-1
                    if sum(LLTh(r-nhood:r+nhood)) == nhood*2+1
                        LLTh2(r)=1;      
                    else, LLTh2(r) = 0;
                    end
                end
                for r = GG-1:GG
                    if sum(LLTh(r-nhood:GG)) == length(LLTh((r-nhood:GG)))
                        LLTh2(r) = 1;         
                    else, LLTh2(r) = 0;
                    end
                end
            
        % Dilate data to restore eroded points - half-window size = 1
            LLTh3 = zeros(size(LLTh));
                for r = 1
                    if sum(LLTh2(1:(r+nhood))) >= 1
                        LLTh3(r) = 1;      
                    else, LLTh3(r) = 0;
                    end
                end
                for r = 2:GG-1
                    if sum(LLTh2(r-nhood:r+nhood)) >=1
                        LLTh3(r)=1;          
                    else, LLTh3(r) = 0;
                    end
                end
                for r = GG-1:GG
                    if sum(LLTh2(r-nhood:GG)) >= 1
                        LLTh3(r) = 1; 
                    else, LLTh3(r) = 0;
                    end
                end
        
        % Dyadic upsampling at odd indices must be performed once for each
        % decomposition level.
            dse0 = LLTh3;
            dse1 = dyadup(dse0,1); %dyadic upsampling, odd indices (1)
            dse2 = dyadup(dse1,1); %dyadic upsampling, odd indices (2)
            dse3 = dyadup(dse2,1); %dyadic upsampling, odd indices (3)
            dse4 = dyadup(dse3,1); %dyadic upsampling, odd indices (4)
          %  dse5 = dyadup(dse4,1); 

           % if DC = 4,
           %     dse_U = dse4;
           % elseif DC = 5
           %     dse_U = dse5;
           % end
            
        % close events - half-window = half * decomp level squared to close
        % the effect of dyadic upsampling
            GG = length(dse4);
            nhood = 0.5*dc; 
            dse4_c = zeros(GG,2);
                for r = 1:nhood
                    if sum(dse4(1:(r+nhood))) ~= 0
                        dse4_c(r) = 1;
                        
                    else, dse4_c(r) = 0;
                    end
                end
                for r = nhood+1:GG-nhood
                    if sum(dse4(r-nhood:r+nhood)) ~=0
                        dse4_c(r)=1;
                    else, dse4_c(r) = 0;
                    end
                end
                for r = GG-nhood+1:GG
                    if sum(dse4(r-nhood:GG)) ~=0
                        dse4_c(r)=1;
                        
                    else, dse4_c(r) = 0;
                    end
                end

        % Erode the ends of events
                dse4_f = zeros(size(dse4_c));
                for r = 1:nhood
                    if sum(dse4_c(1:(r+nhood))) == ...
                            length(dse4_c(1:(r+nhood)))
                        dse4_f(r) = 1;
                    else ,dse4_f(r) = 0;
                    end
                end
                for r = nhood+1:GG-nhood
                    if sum(dse4_c(r-nhood:r+nhood)) == nhood*2+1
                        dse4_f(r)=1;
                    else ,dse4_f(r) = 0;
                    end
                end
                for r = GG-nhood+1:GG
                    if sum(dse4_c(r-nhood:GG)) == ...
                            length(dse4_c((r-nhood):GG))
                        dse4_f(r) = 1;
                    else ,dse4_f(r) = 0;
                    end
                end
        
            LL_above_th = find(dse4_f==1); 
                % pick out line lengths greater than 
                % BaseMedian + SF*StDev. return as 
                % row index (time) of event start
                       
            clear LLAve LLAveSD LLTh LLTh2 LLTh3 nhood GG r dse0 dse1 ...
                dse2 dsde3 dse4 dse4_c 
 %%               
    if isempty(LL_above_th)~=1
        
        %% Event classification    
            % identify hits that are greater than 1 from neighboring event
                evnt_start = zeros(length(LL_above_th),1);  
                    for r = 2:length(LL_above_th)
                        evnt_start(r,1) = LL_above_th(r) - LL_above_th(r-1);
                    end
                evnt_indx = find(evnt_start ~= 1);
       
        % Find event lengths
                evnt_lngth = zeros(length(evnt_indx),1);
                    for r = 1:length(evnt_indx)-1
                        evnt_lngth(r,1) = LL_above_th(evnt_indx(r+1)-1) ...
                            - LL_above_th(evnt_indx(r));
                    end
                evnt_lngth(length(evnt_indx),1) = LL_above_th(end,1)...
                    - LL_above_th(evnt_indx(end,1));

                [row1, col1] = size(evnt_lngth);

        % Find event ends
                evnt_end = zeros(size(evnt_indx));
                    for r = 1:length(evnt_indx-1)
                        evnt_end(r,1) = LL_above_th(evnt_indx(r))+ evnt_lngth(r);
                    end
                ev_end = [LL_above_th(evnt_indx) evnt_end];
              %  ev_end(1,3) = 400;
                ev_end(1,3) = f;
                    for r = 2:length(evnt_indx)
                            ev_end(r,3) = ev_end((r),1)-ev_end((r-1),2);
                    end
    
        % Identify events less than 500 ms from neighbor.
        % Events closer than 500 ms to neighbor will be combined with neighbor.
                ev_end(:,4) = 0;
                    for r=1:length(evnt_indx)
                        if ev_end(r,3) <= .5*f
                           ev_end(r,4) = 0;
                        else ,ev_end(r,4) = 1;
                        end
                    end
                n_to_combine = find(ev_end(:,4) == 1); 
                    % index of events closer than 250 ms
    
        % recreate event list after combining neighbors 
                events = zeros(length(n_to_combine),2);
                    for r = 1:length(n_to_combine)-1
                        events(r,1) = ev_end(n_to_combine(r),1);
                        events(r,2) = ev_end((n_to_combine(r+1)-1),2);
                    end
                events(length(n_to_combine),1) = ...
                    ev_end(n_to_combine(end),1);
                events(length(n_to_combine),2) = ev_end(end,2);
    
        % Sort events based on event lengths
            % Short = events < 5 sec, includes spikes and abnormal events; 
            % Seizures = events > 5 sec in length 
                events(:,3) = events(:,2)-events(:,1);
                [rows, columns] = size(events);
                    for r = 1:rows
                        %if events(r,3) < 5*f
                         if events(r,3) < 5*f
                            events(r,4) = 1;
                        else ,events(r,4) = 2;
                        end
                    end
                ev1 = find(events(:,4)==1);
                ev2 = find(events(:,4)==2);

                    
                    if isempty(ev2) ~= 1
                        seizures = events(ev2,1:2)/f;
                    else seizures = [];
                    end

        % Sort Short events based on event amplitude
            % Spikes = events > 250 uV in amplitude
            % Abnormal = events < 250 uV in amplitude
                if isempty(ev1) ~= 1 
                    short = events(ev1,1:2);
                    [rows, ~] = size(short);
                    SpMin = zeros(rows,2);
                    
                   EEGev = zeros(length(dse4_f),1);
                   EEGev(1:length(EEG)) = EEG(:,col);
                    for r = 1:rows
                        sp = short(r,1):short(r,2);
                        sp = sp';
                        SpMinMax(r,:) = [min(EEGev(sp)) max(EEGev(sp))];
                    end
        
                SpAmp = SpMinMax(:,2)-SpMinMax(:,1);
                    Abn0 = find(SpAmp < 250); %changed to 250 for BONN data
                        abnormal = short(Abn0,1:2)/f;
                    Sp0 = find(SpAmp >= 250); %changed to 250 for BONN data
                        spikes = short(Sp0,1:2)/f;    
                else
                        spikes = [];
                        abnormal = [];
                end

                spikecount = length(spikes);
                abnormalcount = length(abnormal);
                seizurecount = length(seizures);
                    counts = [spikecount seizurecount abnormalcount ];
          
     % find the actual abnormal time totals
    
                if isempty(spikes) ~=1
                    TotSPTime_byAlg = sum(spikes(:,2)-spikes(:,1));
                else ,TotSPTime_byAlg = 0;
                end
                if isempty(abnormal)~=1
                    TotABTime_byAlg = sum(abnormal(:,2)-abnormal(:,1));
                else ,TotABTime_byAlg = 0;
                end
                if isempty(seizures) ~= 1
                    TotSZTime_byAlg = sum(seizures(:,2)-seizures(:,1));
                else ,TotSZTime_byAlg = 0;
                end
                    RealTotTime = [TotSPTime_byAlg  ...
                        TotSZTime_byAlg TotABTime_byAlg];
           %%         
                    clear evnt_start LL_above_th evnt_indx evnt_lngth row1...
                        col1 evnt_end ev_end n_to_combine events ev1 ev2 ...
                        SpMinMax SpAmp Abn0 Sp0 spikecount abnormalcount ...
                        seizurecount TotSPTime_byAlg TotABTime_byAlg ...
                        TotSZTime_byAlg sp short dse4_f dse5
         
        else 
            seizures = [];
            abnormal  = [];
            spikes  = [];
            counts = [0 0 0];
            RealTotTime = [0 0 0];
    end

    
   
    
            % fill in structure with signal analysis
    
                SigSumm(col,k).File = name;
                SigSumm(col,k).TotRecordLength = length(EEG)/f;
               % SigSumm(k).BulkSigScore = Thresholded(k);
                SigSumm(col,k).EventCounts = counts;
                SigSumm(col,k).Spikes = spikes;
                SigSumm(col,k).Seizure = seizures;
                SigSumm(col,k).Abnormal = abnormal; 
                SigSumm(col,k).EventTimeTotals = RealTotTime;
    
                %%
        
     %plot goes here.
         figure
        subplot(2,1,1)
        plot(EEG)
        
              
        t = length(EEG/173.61); 
        eeg_plot = plot((1/173.61:1/173.61:t/173.61),EEG,'k');
        title(name)
        ax = gca;
        ax.Box = 'off';
        ax.YLabel.String = 'µV';
        ax.XLabel.String = 's';
        axis([0 24 -2000 1000])
        axis('on')
        
        hold on
        
        
         bn = 1;
                
                if length(seizures) > 0
                
                    [r c] = size(seizures);
                    
                for b = 1:r
                    fs = seizures(b,1):1/f:seizures(b,2);
                    fl = length(fs);
                    fsz(bn:bn+fl-1) = fs;
                    bn = bn+fl;
                end
                
                sz_plot = plot(fsz,-1000,'r.');
                 
                clear bn fs fl b r c fsz sz_plot
                
                else
                end
        
       
                
                  if length(spikes) > 1
                
                bn = 1;
                [r c] = size (spikes);
                
                for b = 1:r
                    fs = spikes(b,1):1/f:spikes(b,2);
                    fl = length(fs);
                    fsp(bn:bn+fl-1) = fs;
                    bn = bn+fl;
                end
                
                 sp_plot = plot(fsp, -1200, 'c.');
                 
                 clear bn fs fl b fsp sp_plot
                 
                else 
                  end
                
      
                
                     if length(abnormal) > 1
                
                bn = 1;
                [r c] = size (abnormal);
                
                for b = 1:r
                    fs = abnormal(b,1):1/f:abnormal(b,2);
                    fl = length(fs);
                    fabn(bn:bn+fl-1) = fs;
                    bn = bn+fl;
                end
                
                 abn_plot = plot(fabn, -1400, 'g.');
                 
                 clear bn fs fl b fsp abn_plot
                 
                else 
                end
      
        
       
       legend('SF = .5')
       
    %saveas(gcf,['/Users/bergstromr/desktop/Bonn_Dataset/171128/SFFigs/SF005/' name '_SF005.eps'])
     print (gcf,['/Users/bergstromr/desktop/Bonn_Dataset/171128/SFFigs/SF2_5/' name '_SF3_5'],'-depsc2')   
     % saveas(gcf,['/Users/bergstromr/desktop/Bonn_Dataset/171128/SFFigs/SF005/' name '_SF005.pdf'])
    
    close(gcf)
   
       
       %%
                              
        clear spHitssp spHitsabn spHitssz abnHitssp abnHitsabn ...
            abnHitssz szHitssp szHitsabn szHitssz spFNeg abnFNeg ...
            szFNeg SpikeFNeg AbnFNeg SeizureFNeg spFPos abnFPos ...
            szFPos SpikeFNeg AbnFNeg SeizureFNeg SpikeFPos ...
            AbnFPos SeizureFPos FNEGs FPOSs n_n spnc mdnc sznc spc mdc ...
            szc RealVid VidFNeg VidFPos_SP VidFPos_Abn VidFPos_SZ N_N ...
            ConfMat ConfMat2 SMinMax
        
    end

        clear SpEye AbnEye SzEye ByEye

end  
%% 


    csvname =  [name '_SigSumm'];
       % any2csv(SigSumm, ',',1, [csvname '.csv']);
      writetable(struct2table(SigSumm), ...
          ['/Users/bergstromr/desktop/Bonn_Dataset/171128/SFTables/' csvname '_SF3_5.xlsx'])

 
      disp ('donezo')
      
      toc
    