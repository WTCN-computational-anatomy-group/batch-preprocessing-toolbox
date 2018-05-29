function show_img(X,Y,lambda,sd,modality,fn)

C   = numel(Y);
dm1 = size(Y{1});

% Set-up figure
fh = findobj('Type','Figure','Name',fn);
if ~isempty(fh)
    fig = fh;
else
    fig = figure('Name',fn,'NumberTitle','off'); 
end
set(0,'CurrentFigure',fig);  

for c=1:C % Loop over channels
    subplot(C,2,2*(c - 1) + 1)
    if strcmpi(modality,'CT')
        dm = size(X{c});
        imagesc(X{c}(:,:,floor(dm(3)/2) + 1)',[0 100]);
    else
        if iscell(X{c})
            dm = size(X{c}{1});
            imagesc(X{c}{1}(:,:,floor(dm(3)/2) + 1)'); % Show only first within-channel image
        else
            dm = size(X{c});
            imagesc(X{c}(:,:,floor(dm(3)/2) + 1)');
        end
    end
    axis off image; colormap(gray);       
    if iscell(X{c})
        title(['sd=' num2str(round(sd{c}(1),2))])
    else
        title(['sd=' num2str(round(sd(c),2))])
    end

    subplot(C,2,2*(c - 1) + 2)
    if strcmpi(modality,'CT')
        imagesc(Y{c}(:,:,floor(dm1(3)/2) + 1)',[0 100]); 
    else
        imagesc(Y{c}(:,:,floor(dm1(3)/2) + 1)'); 
    end
    axis off image; colormap(gray);    
    title(['lambda=' num2str(round(lambda(c),2))])    
end
drawnow
%==========================================================================  