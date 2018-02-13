%{ 
curvature.m takes .avi video files of larva as input, analyzes the curvature, 
and outputs images with curvature outlined and data tables with larva 
Curvature Index values.

Requires fitcircle.m

Last updated 3/11/17
C. Sam Qian
Grueber Lab
Columbia University

%}

%clear the screen and workspace
clear all;
clc;
warning off;

%ask user to select video files
[file,folder] = uigetfile('*','Select .avi files','Multiselect','on');

%if only one file is selected (if 'file' is a string), convert to cell to
%be consistent with selecting multiple files
if isstr(file)==1;
    file={file};
end
numFiles=size(file,2);

%loop each video file
for n=1:numFiles;
    
    %clear variables before processing the next video
    close all;
    clearvars -except file folder numFiles n
    
    %extract specific file name and load video
    [pathstr,filename,ext]=fileparts(fullfile(folder,file{n}));
    imgObj = VideoReader(fullfile(folder, file{n}));
    numFrames=imgObj.NumberofFrames;
    
    time=zeros(numFrames,1);
    RadiusList=zeros(numFrames,300);

    mkdir('output');
    mkdir('output/',filename);
    mkdir('raw');
    mkdir('raw/',filename);
    mkdir('matrix_data');
    
    %manual override what frame to stop at
    numFrames=100;
    
    %for current video, process the specific range of frames
    for i=1:numFrames;
        
        %--------------------IMAGE PROCESSING----------------------------
%       img=imread(sprintf('_raw_manually_corrected_specks/%s/img%d.tif', filename, i));
        img=read(imgObj, i);
        
        %convert to grayscale
        img=rgb2gray(img);
        imwrite(img,sprintf('raw/%s/img%d.tif', filename, i));
        
%         %if image has been edited, restructure it to be processable
%         if size(img,3)~=1
%             img=img(:,:,1:3);
%         end

        %convert black and white. using 0.1 instead of graythresh'd value works well
        %threshold=graythresh(img);
        img2=im2bw(img,0.1); %baseline 0.1

        %remove all object containing fewer than 10 and more than 5000 pixels
        imgsmallremoved = bwareaopen(img2,200);
        imglargeonly = bwareaopen(img2,7000);
        img2 = imgsmallremoved - imglargeonly;

        %to determine size of larva in pixels
        %Total_White_Pixels = nnz(img2); % total area based on number of white pixels
        
        % fill a gap in the pen's cap (10 pixels)
        se = strel('disk',10); %baseline 10
        img2 = imclose(img2,se);

        % fill any holes, so that regionprops can be used to estimate
        % the area enclosed by each of the boundaries
        img2 = imfill(img2,'holes');

        %second round of filtering by size (smaller than 300 pixels and
        %bigger than 1500 pixels)
        imgsmallremoved = bwareaopen(img2,600);
        imglargeonly = bwareaopen(img2,7000);
        img2 = imgsmallremoved - imglargeonly;

        %clear border touching objects
        img2=imclearborder(img2);

        %save bw image in subdirectory
        %imwrite(~img2,sprintf('bw/%s/img%d.tif', filename, i));
        
        %get boundaries
        [B,L] = bwboundaries(img2,'noholes');
        
        %if there are boundary values
        if size(B,1)~=0;
            B_vals=B{1,1};

            %resize boundaries to 300 boundary points
            B_vals=imresize(B_vals,[300 2]);

            %create temporary new array where 50 points are tacked on to each
            %ends of the 300 points array (so the array is circular for angle
            %calculations)
            C_vals(1:50,1:2)=B_vals(251:300,1:2);
            C_vals(51:350,1:2)=B_vals(1:300,1:2);
            C_vals(351:400,1:2)=B_vals(1:50,1:2);

            %Put tripoint data into 3 dimensional array (tripoints, x and y values, edge point)
            c_factor=10;
            for p=1:300;
                Tri(1,1,p)=C_vals(p+50,2);
                Tri(1,2,p)=C_vals(p+50,1);
                Tri(2,1,p)=C_vals(p+50-c_factor,2);
                Tri(2,2,p)=C_vals(p+50-c_factor,1);
                Tri(3,1,p)=C_vals(p+50+c_factor,2);
                Tri(3,2,p)=C_vals(p+50+c_factor,1);

                %midpoint coordinate
                Mid_pt(p,1:2)=[(Tri(2,1,p)+Tri(3,1,p))/2, (Tri(2,2,p)+Tri(3,2,p))/2];
                Mid_pt(p,1:2)=round(Mid_pt(p,1:2),0);

                %if midpoint lies within larva, set column 3 of midpoint to -1,
                %if outside larva, set column 3 to 1
                if img2(Mid_pt(p,2),Mid_pt(p,1))==1
                    Mid_pt(p,3)=-1;
                else
                    Mid_pt(p,3)=1;
                end

                %get circle fit and record inversed radius (w/ inside outside
                %sign multiplier 
                [z, r] = fitcircle(Tri(:,:,p)');
                C_radius(p)=Mid_pt(p,3)/r;
                
                %Optional visualization of curvature analysis
                %{
                imshow(img2);
                hold on;
                %plot circle
                th = 0:pi/50:2*pi;
                xunit = r * cos(th) + z(1);
                yunit = r * sin(th) + z(2);
                plot(xunit, yunit,'c');       

                %add tripoints
                plot(Tri(:,1,p),Tri(:,2,p),'c*')
                
                %add midpoint (green if inside, red if outside animal)
                if Mid_pt(p,3)==1
                    plot(Mid_pt(p,1),Mid_pt(p,2),'g*');
                else
                    plot(Mid_pt(p,1),Mid_pt(p,2),'r*');
                end
                
                text(50,50,['Radius:' num2str(C_radius(p))],'Color', 'w');
                
                pause(0.1);
                hold off;
                %}

            end

                %plot edge points
                figure('visible','off');
                test=repmat(img, [1, 1, 3]);

                %restruct if necessary
                if size(test,3)~=1
                    test=test(:,:,1:3);
                end

                imshow(test);
                colormap(gray);
                hold on;
                scatter(Tri(1,1,:),Tri(1,2,:), 20, C_radius(:), 'filled');
                colormap(jet);
                set(gca,'CLim',[-0.15 0.15]);

                set(gcf, 'InvertHardCopy', 'off');
                saveas(gcf,sprintf('output/%s/img%d.tif', filename, i));

                hold off;
                close(gcf);

                string = sprintf('%s analyzed frame:%d', filename, i);
                disp(string)

                %save the inversed radius (aka Curvature Index)
                RadiusList(i,1:300)=C_radius;
        end
    end
    
    RadiusList(isnan(RadiusList))=0;
    save(sprintf('matrix_data/%s_radius.mat', filename), 'RadiusList');
end;

sprintf('Completed')

