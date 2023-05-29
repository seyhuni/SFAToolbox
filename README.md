# SFAToolbox

## Main folders: :file_folder: 

*GUI:*
* SFAToolbox.m :  
* SFAToolbox.fig : 

* Filtering: 

```
Band    = (2 / SR) * [45, 55];
[B, A]  = butter(4, Band, 'stop');
block_duration_datacorrected_last = filtfilt(B, A, double(block_duration_datacorrected));
N=length(block_duration_datacorrected_last);         %number of points
t=(0:N-1)/SR;   %time vector
sgf = sgolayfilt(block_duration_datacorrected_last,3,201);
ynew=sgf/max(sgf); 
% initialize filtered signal
eogF = ynew;
```

Teager–Kaiser energy operator to obtain EMG Bursts:
```
for i=2:length(eogF)-1
    eogF(i) = ynew(i)^2 - ynew(i-1)*ynew(i+1);
end   
    
```

Wavelet transform to get swallow patterns and RMS for burst detection:
```
[c,l] = wavedec(eogF,7,'sym4'); % 8 level decomposition  
for t=1:7
    A(:,t)=wrcoef('a',c,l,'sym4',t);
end
A8=A(:,7);
rmsSwallows = sqrt(movmean(A8_.^2, 2500));
```
