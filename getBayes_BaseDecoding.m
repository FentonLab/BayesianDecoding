%compute bayesian decoding and compare errors before-after PCP
%setup decoding from base
clear; close all; clc;

addpath('/Volumes/DATA/matlab/my_functions/');

minSpikes = 100;
minCells = 5;
eegFS = 2000;

NedgesXY = 32;
%edgesXY = 0.5:8:256.5;
%NedgesXY = length(edgesXY)-1;

winLenHalfSpk = 150;


IND = 1:NedgesXY*NedgesXY;
[iY,iX] = ind2sub([NedgesXY,NedgesXY],IND);

drTetrodeTS = 'TETRODE_TS/';

drART = 'ART_border025dif/';

drXY = 'XY/';

%cell lists, file lists
load('code phaseLFP/filesCellsPC.mat');
%load('code phaseLFP/filesCells.mat');

%plotDose = 3; %0=SAL,3=3mg,5=5mg
doses = [0,3,5];

loadData = 0;

coefXY = 82/32;

colors = [0.3 0.3 0.3; 0 0 1; 1 0 0];

if loadData == 1

    resultsErr = cell(1,3);
    Ncells = cell(1,3);
    Posteriors = cell(1,3);
    Errors = cell(1,3);
    Locations = cell(1,3);

    for doseI = 1:3

        plotDose = doses(doseI);

        errs = nan(length(files),2);
        Nc = nan(length(files),2);
        postDose = cell(length(files),2);
        errorsDose = cell(length(files),2);
        locDose = cell(length(files),2);

        for gI = 1:length(files)

            fn = files{gI,2};

            %get cellList 
            pairsBase1 = filesCells{gI,1}; %cells from BASE (first file)
            pairsBase2 = filesCells{gI,2}; %cells from BASE (first file)

            if size(pairsBase1,1) ~= size(pairsBase2,1); disp(['no match' fn]); continue; end;
            d = pairsBase1-pairsBase2;
            if sum(d(:)) ~= 0; disp('different'); continue; end;

            %cells
            pairsBase = filesCells{gI,1}; %cells from BASE (first file)
            hBase = cell(size(pairsBase,1),2,2); %need to load for each file
            mouseID = [];
            expIDs = cell(1,3);

            %get type
            if exist([drTetrodeTS fn],'file') == 0; continue; end;
            s = strsplit(fn,'-');
            s = s{2};

            if plotDose == 0 && strcmp(s,'SAL') ~= 1; continue; end;
            if plotDose == 3 && strcmp(s,'PCP3') ~= 1; continue; end;
            if plotDose == 5 && strcmp(s,'PCP5') ~= 1; continue; end;

            Px = [];
            FRMS = [];

            for fI = 1:2

                fn = files{gI,fI};
                if exist([drTetrodeTS fn],'file') == 0; continue; end;
                disp(fn);
                load([drTetrodeTS fn]);
                load([drXY fn]);

                %convert time to ms
                roomTimeStamps = roomTimeStamps';
                tetrodeTimestamp = double(tetrodeTimestamp)/10;
                roomTimeStamps = double(roomTimeStamps)/10;

                %roomTimeStamp might have steps
                k = find(diff(roomTimeStamps) < 0);
                for i = 1:length(k)
                    roomTimeStamps(k(i)+1:end) = roomTimeStamps(k(i)+1:end) + roomTimeStamps(k(i));
                end

                s = strsplit(fn,'-');
                expIDs{fI} = s{2};
                if fI == 1
                    mouseID = s{1};
                end

                %find X,Y for all tetrodeTS
                edges = mean([roomTimeStamps(2:end); roomTimeStamps(1:end-1)]);
                d = edges(10)-edges(9);
                edges = cat(2,edges(1)-d, edges, edges(end)+d);
                ind = discretize(tetrodeTimestamp,edges);

                knan = ~isnan(ind);
                tetrodeTimestamp = tetrodeTimestamp(knan);
                tetrodeChannel = tetrodeChannel(knan);
                tetrodeUnit = tetrodeUnit(knan);
                ind = ind(knan);

                X = positionXY(1,:);
                Y = positionXY(2,:);

                kbad = X == 0 & Y == 0;
                X(kbad) = nan;
                Y(kbad) = nan;

                Xspk = X(ind);
                Yspk = Y(ind);

                %check min spikes
                kKeep = true(1,size(pairsBase,1));
                for cellI = 1:size(pairsBase,1)
                    T = pairsBase(cellI,1);
                    U = pairsBase(cellI,2);
                    k = tetrodeChannel == T & tetrodeUnit == U;
                    if sum(k) < minSpikes; kKeep(cellI) = 0; end;
                end
                pairsBase = pairsBase(kKeep,:);

                if size(pairsBase,1) < minCells; disp(['few cells' fn ' ' num2str(size(pairsBase,1))]); continue; end;

                Nc(gI,fI) = size(pairsBase,1);

                Nwin = length(roomTimeStamps);

                %minmax values of XY
                minX = min(X);
                maxX = max(X);
                minY = min(Y);
                maxY = max(Y);

                %edgesXY = 0.5:8:256.5;
                edgesX = linspace(minX-0.1,maxX+0.1,NedgesXY+1);
                edgesY = linspace(minY-0.1,maxY+0.1,NedgesXY+1);

                %dwell distribution
                %spatial occupancy
                kGoodSpot = ~isnan(X) & ~isnan(Y);
                nnXY = histc3(cat(2,Y(kGoodSpot)',X(kGoodSpot)'),'Edges',{edgesY,edgesX});
                nnXY = nnXY(1:end-1,1:end-1);
                %convert to time
                nnXYt = nnXY*0.033; %33ms between frames

                %probability of place
                if fI == 1
                    Px = nnXY / nansum(nnXY(:));
                end

                %bin position
                XYS = nan(length(X),2);
                [~,indX] = histc(X,edgesX);
                [~,indY] = histc(Y,edgesY);
                XYS(:,1) = indX;
                XYS(:,2) = indY;

                XYS(~kGoodSpot,:) = nan;

                locDose{gI,fI} = XYS;


                if fI == 1
                    %firing rate maps
                    FRMS = cell(1,size(pairsBase,1));

                    %cells
                    for cellI = 1:size(pairsBase,1)

                        t = pairsBase(cellI,1);
                        u = pairsBase(cellI,2);

                        k = tetrodeChannel == t & tetrodeUnit == u;
                        if sum(k) == 0; disp('no spikes!'); return; end;

                        %Firing rate map
                        Xu = Xspk(k);
                        Yu = Yspk(k);

                        %number of spikes at given pixel
                        nn = histc3(cat(2,Yu',Xu'),'Edges',{edgesY,edgesX});
                        nn = nn(1:end-1,1:end-1);

                        %firing rate map
                        frm = nn./nnXYt;

                        %smooth with gaussian kernel, size, std
                        h = fspecial('gaussian', [3,3], 1);
                        frmS = filter2(h,frm);

                        FRMS{cellI} = frmS;
                    end
                end

                %firing rates
                FRS = zeros(size(pairsBase,1),Nwin);

                %cells
                for cellI = 1:size(pairsBase,1)

                    t = pairsBase(cellI,1);
                    u = pairsBase(cellI,2);

                    k = tetrodeChannel == t & tetrodeUnit == u;
                    if sum(k) == 0; disp('no spikes!'); return; end;

                    tsu = tetrodeTimestamp(k);

                    %Firing rate
                    for wI = 1:Nwin
                        st = roomTimeStamps(wI)-winLenHalfSpk;
                        ed = roomTimeStamps(wI)+winLenHalfSpk;
                        FRS(cellI,wI) = sum(tsu>st & tsu<ed);
                    end

                    %Firing rate map
                    Xu = Xspk(k);
                    Yu = Yspk(k);

                    %number of spikes at given pixel
                    nn = histc3(cat(2,Yu',Xu'),'Edges',{edgesY,edgesX});
                    nn = nn(1:end-1,1:end-1);

                    %firing rate map
                    frm = nn./nnXYt;

                    %smooth with gaussian kernel, size, std
                    h = fspecial('gaussian', [3,3], 1);
                    frmS = filter2(h,frm);

                    FRMS{cellI} = frmS;
                end

                FRSsum = sum(FRS,1);


                %compute Bayes from maxima
                %Bayes  
                X1 = Px;
                XYp = nan(Nwin,2); %predicted
                Nspk = zeros(Nwin,1);
                Pxns = nan(NedgesXY,NedgesXY,Nwin);
                for wI = 1:Nwin
                    %cells
                    X2 = ones(NedgesXY);
                    X3 = zeros(NedgesXY);
                    for uI = 1:size(pairsBase,1)
                        ni = FRS(uI,wI); %number of spikes
                        Nspk(wI) = Nspk(wI) + (ni > 0);
                        X2 = X2 .* (FRMS{uI}.^ni);
                        X3 = X3 + FRMS{uI};
                    end
                    X3 = X3 * (-2*winLenHalfSpk/1000); %delta in seconds
                    X3 = exp(1).^X3;

                    Pxn = X1.*X2.*X3;

                    %maximum
                    %[maxA,ind] = max(Pxn(:));
                    %[yi,xi] = ind2sub(size(Pxn),ind);

                    %{
                    %prob-weighted average (less noisy than maxima)
                    Pxn = Pxn(:)'; %linearize
                    Pxn = Pxn/nansum(Pxn); %prob
                    xi = nansum(iX.*Pxn);
                    yi = nansum(iY.*Pxn);
                    %}


                    %prob-weighted average (less noisy than maxima)
                    %only for top 10% of pixels
                    %Pxn = Pxn(:)'; %linearize
                    %Pxn = Pxn/nansum(Pxn); %prob

                    [PxnSort,ksort] = sort(Pxn(:),'descend');
                    N = floor(length(Pxn(:))/10);
                    iXsort = iX(ksort);
                    iYsort = iY(ksort);

                    knan = ~isnan(PxnSort); %sorted values start with nan!
                    iXsort = iXsort(knan);
                    iYsort = iYsort(knan);
                    PxnSort = PxnSort(knan);

                    iXsort = iXsort(1:N); %select 10%
                    iYsort = iYsort(1:N);
                    PxnSort = PxnSort(1:N);

                    s = nansum(PxnSort); %sum of probs
                    xi = nansum(iXsort.*PxnSort')/s;
                    yi = nansum(iYsort.*PxnSort')/s;

                    XYp(wI,1) = xi;
                    XYp(wI,2) = yi;

                    %save posterior
                    Pxn = Pxn/nansum(Pxn(:));
                    Pxns(:,:,wI) = Pxn;
                end
                err = ((XYS(:,1) - XYp(:,1)).^2 + (XYS(:,2) - XYp(:,2)).^2).^0.5; %real-predicted

                %err = XYS(:,1) - XYp(:,1);
                errorsDose{gI,fI} = err;

                kNoSpk = (FRSsum == 0);

                err(kNoSpk) = nan;

                k = find(kNoSpk);
                for i = 1:length(k)
                    ki = k(i);
                    Pxns(:,:,ki) = nan;
                end

                postDose{gI,fI} = Pxns;


                %figure; hold on;
                %plot(XYS(~kNoSpk,1),'k');
                %plot(XYp(~kNoSpk,1),'r');

                errs(gI,fI) = nanmean(err);

            end %files

        end %groups

        %average error difference
        e = errs(:,1) - errs(:,2);

        k = ~isnan(e);
        e = e(k);

        resultsErr{doseI} = e;

        Ncells{doseI} = Nc(k,1);

        Posteriors{doseI} = postDose;

        Errors{doseI} = errorsDose;

        Locations{doseI} = locDose;

    end

    save('bayes.mat','resultsErr','Ncells','Posteriors','Errors','Locations','-v7.3');
else
    load('bayesPlace150.mat');
end

return;

%export stats for error average
labels = {'PCP0','PCP3','PCP5'};
fid = fopen('stats.dat', 'w');
fprintf(fid, '%s\t%s\n', 'dose','err');
for i = 1:3
    x = resultsErr{i};
    x = x(~isnan(x));
    for j = 1:length(x)
        fprintf(fid, '%s\t%12.8f\n', labels{i},x(j));
    end
end
fclose(fid);

%show number of cells
for doseI = 1:3
    x = Ncells{doseI}
end
    

%test normality for error from all time windows
for i = 1:3
    disp(['======== ' num2str(i) ' ========']); 
    x = Errors{i};
    for j = 1:size(x,1)
        if isempty(x{j,2}); continue; end;
        xj = x{j,2};
        [H,P,KSSTAT,CV] = kstest(xj);
        disp(num2str(P));
    end
end

%test normality for error averages
for i = 1:3
    disp(['======== ' num2str(i) ' ========']); 
    x = resultsErr{i};
    [H,P,KSSTAT,CV] = kstest(x);
    disp(P);
end

%plot error averages
figure('Position',[100 100 400 300]); hold on;
for doseI = 1:3
    x = resultsErr{doseI};
    x = x(~isnan(x));
    x = x*coefXY; %make result in cm
    x = -1*x; %reverse
    m = nanmean(x);
    s = nanstd(x)/sqrt(length(x));
    bar(doseI, m, 'FaceColor',colors(doseI,:),'EdgeColor','none');
    errorbard(doseI,m,s,2,0.2,colors(doseI,:))
end
ylabel('Error difference [cm]'); 


%average error size
%82cm = 32pix, 1pix = 2.56cm
coef = 2.56;
for i = 1:3
    disp(['======== ' num2str(i) ' ========']); 
    x = Errors{i};
    errAver = [];
    for j = 1:size(x,1)
        if isempty(x{j,2}); continue; end;
        xj = x{j,2};
        xj = nanmean(xj) * coef;
        errAver = cat(1,errAver,xj);
    end
    disp(errAver)
end

%average error size
%82cm = 32pix, 1pix = 2.56cm
coef = 2.56;
for i = 1:3
    disp(['======== ' num2str(i) ' ========']); 
    x = Errors{i};
    errAver = [];
    for j = 1:size(x,1)
        if isempty(x{j,2}); continue; end;
        xj = x{j,2};
        xj = nanmean(xj) * coef;
        errAver = cat(1,errAver,xj);
    end
    disp(nanmean(errAver))
end

%centered posterior
labels = {'PCP 0','PCP3','PCP5'};
pAvers = cell(1,3);
for doseI = 1:3
    x = Posteriors{doseI}; %posteriors
    l = Locations{doseI}; %locations
    pAver = zeros(66,66); %final average
    pCenter = 33; %middle
    Naver = 0; %number of averaged posteriors
    for expI = 1:size(x,1) %experiments
        xExp = x{expI,2}; %post
        if isempty(xExp); continue; end;
        if sum(~isnan(xExp(:))) == 0; continue; end;
        lExp = l{expI,2}; %location
        %posteriors
        for nI = 1:size(xExp,3)
            pxx = xExp(:,:,nI); %post
            if sum(~isnan(pxx(:))) == 0; continue; end; %ignore if empty
            pxx(isnan(pxx)) = 0; %convert nans to zeros
            loc = lExp(nI,:); %location
            X = loc(1);
            Y = loc(2);
            if isnan(X) || isnan(Y); continue; end;
            Xp = pCenter - X;
            Yp = pCenter - Y;
            pAver(Yp:Yp+32-1,Xp:Xp+32-1) = pAver(Yp:Yp+32-1,Xp:Xp+32-1) + pxx;
            Naver = Naver + 1;
        end
    end
    pAver = pAver / Naver; %average
    
    pAvers{doseI} = pAver;
end
figure('Position',[100 100 1000 300]); hold on;
for doseI = 1:3
    subplot(1,3,doseI); 
    %no log
    x = pAvers{doseI};
    imagesc(x,[0 0.001]);
    %log
    %x = log(pAvers{doseI});
    %imagesc(x,[-10 -6]);
    colormap(parula);
    axis ij tight equal
    title(labels{doseI});
    colorbar
end
    

%size of posterior
figure('Position',[100 100 400 300]); hold on;
%export stats for posterior size
labels = {'PCP0','PCP3','PCP5'};
fid = fopen('stats_posterior_size_150x.dat', 'w');
fprintf(fid, '%s\t%s\n', 'dose','sizeRatio');
postSizes = cell(3,6,2);
for doseI = 1:3
    x = Posteriors{doseI};
    
    diffs = [];
    
    counter = 1;
    
    %exp
    for expI = 1:size(x,1)
        xB = x{expI,1}; %before injection
        xA = x{expI,2}; %after
        
        if isempty(xB); continue; end;
        if isempty(xA); continue; end;
        
        %xB(xB == 0) = nan;
        %xA(xA == 0) = nan;
        
        %xB = log(xB);
        %xA = log(xA);

        %get average pixel values for threshold - before only
        xBa = xB(:);
        m = nanmean(xBa);
        s = nanstd(xBa);
        th = m+s;
        
        %get numbers of pixels above threshold
        nB = nan(1,size(xB,3));
        nA = nan(1,size(xA,3));
        for nI = 1:size(xB,3)
            pxx = xB(:,:,nI);
            pxx = pxx(:);
            nB(nI) = sum(pxx>th);
        end
        for nI = 1:size(xA,3)
            pxx = xA(:,:,nI);
            pxx = pxx(:);
            nA(nI) = sum(pxx>th);
        end
    
        postSizes{doseI,counter,1} = nB; %store posterior sizes
        postSizes{doseI,counter,2} = nA;
        counter = counter + 1;
        
        %difference of after/before
        d = nanmean(nA)/nanmean(nB);
        
        diffs = cat(1,diffs,d);
    end
    
    for i = 1:length(diffs)
        fprintf(fid, '%s\t%12.8f\n', labels{doseI},diffs(i));
    end
    
    m = nanmean(diffs);
    s = nanstd(diffs)/sqrt(length(diffs));
    bar(doseI, m, 'FaceColor',colors(doseI,:),'EdgeColor','none');
    errorbard(doseI,m,s,2,0.2,colors(doseI,:))
end
fclose(fid);

%plot distributions of posterior sizes
figure('Position',[100 100 700 300]); hold on;
col = [0.5 0.5 0.5; 1 0 0];
edges = 0:10:50;
for doseI = 1:3
    
    %experiments
    for expI = 1:6
        
        subplot(3,6,(doseI-1)*6 + expI); hold on;
        
        %before after
        for baI = 1:2
            xExp = postSizes{doseI,expI,baI};
            if isempty(xExp); continue; end;
            
            return;
            
            h = histc(xExp,edges);
            h = h/nansum(h);
            plot(edges(1:end-1),h(1:end-1),'Color',col(baI,:),'LineWidth',2);
            axis([edges(1) edges(end) 0 0.4]);
        end
    end
end


%plot distributions of error
figure('Position',[100 100 700 300]); hold on;
col = [0.5 0.5 0.5; 1 0 0];
edges = 0:2:50;
for doseI = 1:3
    x = Errors{doseI}; 
    
    k = ~cellfun(@isempty,x(:,1)) & ~cellfun(@isempty,x(:,2));
    x = x(k,:);
    
    %experiments
    for expI = 1:size(x,1)
        
        subplot(3,6,(doseI-1)*6 + expI); hold on;
        
        %before after
        for baI = 1:2
            xExp = x{expI,baI};
            xExp = xExp*coefXY;
            
            h = histc(xExp,edges);
            h = h/nansum(h);
            plot(edges(1:end-1),h(1:end-1),'Color',col(baI,:),'LineWidth',2);
            axis([edges(1) edges(end) 0 0.11]);
        end
        
    end
end

%plot posterior average
figure('Position',[100 100 400 300]); hold on;
%export stats for posterior size
labels = {'PCP0','PCP3','PCP5'};
%fid = fopen('stats_posterior_size_150.dat', 'w');
%fprintf(fid, '%s\t%s\n', 'dose','sizeRatio');
postVar = cell(3,6,2);
for doseI = 1:3
    x = Posteriors{doseI};
    
    diffs = [];
    
    counter = 1;
    
    %exp
    for expI = 1:size(x,1)
        xB = x{expI,1}; %before injection
        xA = x{expI,2}; %after
        
        if isempty(xB); continue; end;
        if isempty(xA); continue; end;
        
        %compute variances for each posterior
        vB = nan(1,size(xB,3));
        vA = nan(1,size(xA,3));
        for i = 1:size(xB,3)
            xi = xB(:,:,i); xi = xi(:);
            vB(i) = nanmean(xi);
        end
        for i = 1:size(xA,3)
            xi = xA(:,:,i); xi = xi(:);
            vA(i) = nanmean(xi);
        end
        
        %postSizes{doseI,counter,1} = nB; %store posterior sizes
        %postSizes{doseI,counter,2} = nA;
        counter = counter + 1;
        
        %difference of after/before
        d = nanmean(vA)/nanmean(vB);
        
        diffs = cat(1,diffs,d);
    end
    
    for i = 1:length(diffs)
    %    fprintf(fid, '%s\t%12.8f\n', labels{doseI},diffs(i));
    end
    
    m = nanmean(diffs);
    s = nanstd(diffs)/sqrt(length(diffs));
    bar(doseI, m, 'FaceColor',colors(doseI,:),'EdgeColor','none');
    errorbard(doseI,m,s,2,0.2,colors(doseI,:))
end
%fclose(fid);

















%{
%compute error from all pixels
%Bayes  
X1 = Px;
[yi,xi] = ind2sub(size(X1),ind);
err = nan(1,Nwin);
Nspk = zeros(Nwin,1);
for wI = 1:Nwin
    %cells
    X2 = ones(NedgesXY);
    X3 = zeros(NedgesXY);
    for uI = 1:length(USel)
        ni = FRS(uI,wI); %number of spikes
        Nspk(wI) = Nspk(wI) + (ni > 0);
        X2 = X2 .* (FRMS{uI}.^ni);
        X3 = X3 + FRMS{uI};
    end
    X3 = X3 * (-2*winLenHalfSpk/1000); %delta in seconds
    X3 = exp(1).^X3;

    Pxn = X1.*X2.*X3;

    Pxn = Pxn/nansum(Pxn(:)); %normalize

    Pxn = Pxn(:)'; %linearize

    Xo = XYS(wI,1);
    Yo = XYS(wI,2);
    if Xo == 0 || Yo == 0 || isnan(Xo) || isnan(Yo); continue; end;

    eXY = ((Xo - iX).^2 + (Yo - iY).^2).^0.5; %real-predicted

    e = eXY.*Pxn;

    err(wI) = nanmean(e);
end
%}

%{
%remove estimates with no spikes
k = Nspk == 0;
err(k) = nan;
XYp(k,:) = nan;

%err = smth(err,31);

%err = log(err);
%merr = nanmean(err);
%serr = nanstd(err);
%err = (err-merr)/serr;

%}
