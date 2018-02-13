%{ 
kymograph.m takes the matrix files generated from curvature.m as input, 
outputs frame-aligned kymographs with head and tail regions masked with 
negative C.I. values (which are excluded from analysis of positive C.I.s), 
and saves the masked matrix files for quantification.

Last updated 3/11/17
C. Sam Qian
Grueber Lab
Columbia University

%}

%clear the screen and workspace
clear all;
clc;

mkdir('output');
mkdir('matrix_data_masked');

%ask user to select .avi files
[file,folder] = uigetfile('*','Select .avi files','Multiselect','on');

%if only one file is selected (if 'file' is a string), convert to cell to
%be consistent with selecting multiple files
if isstr(file)==1;
    file={file};
end

numFiles=size(file,2);

%loop each file
for n=1:numFiles;
    
    clearvars -except file folder numFiles n
    
    %extract specific file name
    [pathstr,filename,ext]=fileparts(fullfile(folder,file{n}));
    
    %open the matrix file for this video
    load(sprintf('matrix_data/%s_radius.mat', filename));
    
    %get frame number of files in folder
    [presentFrames]=dir(sprintf('raw/%s/*.tif', filename ));
    presentFramesList={presentFrames.name}';
    
    for t=1:size(presentFramesList,1);   
        y=char(presentFramesList(t));
        stringNum(t)=sscanf(y,'img%f.tif',[1 Inf]);
    end
    
    %get the lowest and highest frame number that exists in the raw img folder
    minF=min(stringNum);
    maxF=max(stringNum);
    
    start=minF;
    fin=maxF;

    %% ------- Generating kymographs --------------
    
    %Get time window of roll
    RadiusList=(RadiusList)';
    RadiusList=RadiusList(:,start:fin);
    
    % -------- ALIGNMENT --------------------------
    matrix=RadiusList;
    %repeat from 2nd column to the number of columns
    for j=2:size(matrix,2);
        %shift down by 1 (counter=i) and calculate distance between current and previous column 
        for i=1:size(matrix,1)
            %shift down current column
            matrix(:,j)=circshift(matrix(:,j),1);
            %current column
            A=matrix(:,j);
            %previous column
            B=matrix(:,j-1);
            %get distance 
            D(i,j)=sqrt(sum((A - B) .^ 2));
        end

        %get min distance for current column, j
        [ x, y ] = min(D(:,j));

        %shift current column down by y
        matrix(:,j)=circshift(matrix(:,j),y);
        
%         %at the j'th column, shift by 150 (FOR MANUAL ALIGNMENT ADJUSTMENT)
%          if j==72
%              matrix(:,j)=circshift(matrix(:,j),150);
%          end;
         
        
    end;
    
    % -----------------MASKING----------------------------------
    matrixExpanded(1:50,:)=matrix(251:300,:);
    matrixExpanded(51:350,:)=matrix;
    matrixExpanded(351:400,:)=matrix(1:50,:);

    %for each column
    for a=1:size(matrix,2);
        %for each row
        for b=1:size(matrix,1);
            %get blueness at row i as 10 pixel wide strip
            blueness(b,1)=sum(matrixExpanded(b+50-5:b+50+5,a));
        end

        %get bluest row (blueLoc1)
        [blue1,blueLoc1] = min(blueness);

        %determine how many rows down to circshift to get bluest row to row 150
        rowDif1=150-blueLoc1;

        %In blueness vector, shift bluest row to 150 and zero surrounding
        %regions, then shift it back
        blueness=circshift(blueness,rowDif1);
        blueness(75:225)=0;
        blueness=circshift(blueness,-rowDif1);

        %get second bluest row (blueLoc2)
        [blue2,blueLoc2] = min(blueness);

        rowDif2=150-blueLoc2;

        %--------------- MASKING --------------------------

        maskSize=25;

        %mask region flanking bluest row
        matrix(1:300,a)=circshift(matrix(:,a),rowDif1);
         matrix(150-maskSize:150+maskSize,a)=-0.4;
        
%         uncomment to highlight line only
%         matrix(150-maskSize-1:150-maskSize+1,a)=-1;
%         matrix(150+maskSize-1:150+maskSize+1,a)=-1;
        
        matrix(1:300,a)=circshift(matrix(:,a),-rowDif1);

        %mask region flanking second bluest row
        matrix(1:300,a)=circshift(matrix(:,a),rowDif2);
         matrix(150-maskSize:150+maskSize,a)=-0.4;
        
%         uncomment to highlight line only
%         matrix(150-maskSize-1:150-maskSize+1,a)=-1;
%         matrix(150+maskSize-1:150+maskSize+1,a)=-1;
        
        matrix(1:300,a)=circshift(matrix(:,a),-rowDif2);
    end
    
    %saves updated matrix file
    RadiusList=matrix;
    save(sprintf('matrix_data_masked/%s_radius.mat', filename), 'RadiusList');
    
    %generate kymograph
    imagesc(RadiusList);
    set(gca,'CLim',[-0.15 0.15]);
    set(gcf,'color','w');
    colormap(jet);
    %h=colorbar;
    
    ylabel('Boundary Point position');
    xlabel('Frames');

    %xlabels = get(gca, 'XTickLabel');
    %xlabels = linspace(0,2.5,length(xlabels));
    %set(gca,'XTickLabel',xlabels);

    %save kymograph in output folder
    saveas(gcf,sprintf('output/%s.tif', filename));
    %close(gcf);

end;

sprintf('Completed')
