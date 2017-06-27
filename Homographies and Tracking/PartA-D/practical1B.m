function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

%close all open figures
close all;
tic
%load in the required data
load('PracticalDataSm','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3

%show images and points
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
plot(pts1(1,:),pts1(2,:),'r.'); 
plot(pts1b(1,:),pts1b(2,:),'m.');
figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
plot(pts2(1,:),pts2(2,:),'r.'); 
figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
plot(pts3(1,:),pts3(2,:),'m.'); 

%****TO DO**** 
%calculate homography from pts1 to pts2
H12  = calcBestHomography(pts1,pts2);
H13  = calcBestHomography(pts1b,pts3);
%****TO DO**** 
%for every pixel in image 1
    %transform this pixel position with your homography to find where it 
    %is in the coordinates of image 2
    %if it the transformed position is within the boundary of image 2 then 
        %copy pixel colour from image 2 pixel to current position in image 1 
        %draw new image1 (use drawnow to force it to draw)
    %end
%end;
[ImHeight, ImWidth, ImZ] = size(im1);
[ImHeight2, ImWidth2, ImZ2] = size(im2);
for i = 1:ImWidth
    for j = 1:ImHeight
        homPts2 = H12*[i,j,1]';
        carts = homPts2(1:2, :)/homPts2(3,:);
        carts = round(carts);
        x = carts(1);
        y = carts(2);
        if (all(carts>1) && (x < ImWidth2) && (y< ImHeight2) )
            im1(j,i,:) = im2(y,x,:);
        end
          %drawnow
      
    end
end
%size(carts);

%figure;
%imshow(uint8(im1));
%****TO DO****
%repeat the above process mapping image 3 to image 1.
[ImHeight, ImWidth, ImZ] = size(im1);
[ImHeight2, ImWidth2, ImZ2] = size(im3);
for i = 1:ImWidth
    for j = 1:ImHeight
        homPts3 = H13*[i,j,1]';
        carts = homPts3(1:2, :)/homPts3(3,:);
        carts = round(carts);
        x = carts(1);
        y = carts(2);
        if (all(carts>1) && (x < ImWidth2) && (y< ImHeight2) )
            im1(j,i,:) = im3(y,x,:);
        end
          %drawnow
      
    end
end

figure;
imshow(uint8(im1));
saveas(gcf, 'practical1B', 'jpg')
toc
function H = calcBestHomography(pts1Cart, pts2Cart)

%should apply direct linear transform (DLT) algorithm to calculate best
%homography that maps the points in pts1Cart to their corresonding match in 
%pts2Cart

%****TO DO ****: replace this
H = eye(3);
pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))];
pts2Hom = [pts2Cart; ones(1,size(pts1Cart,2))];
A = [];

for i=1:5
    A = [A ;0, 0, 0, -pts1Hom(:,i)', pts2Hom(2,i)*(pts1Hom(:,i))'];
    A = [A ; pts1Hom(:,i)',0, 0, 0, -pts2Hom(1,i)*(pts1Hom(:,i))'];
    
end
%**** TO DO ****;
%first turn points to homogeneous
%then construct A matrix which should be (10 x 9) in size
%solve Ah = 0 by calling
h = solveAXEqualsZero(A); %(you have to write this routine too - see below)
H = reshape(h,[3,3])';

%reshape h into the matrix H

%Beware - when you reshape the (9x1) vector x to the (3x3) shape of a homography, you must make
%sure that it is reshaped with the values going first into the rows.  This
%is not the way that the matlab command reshape works - it goes columns
%first.  In order to resolve this, you can reshape and then take the
%transpose


%==========================================================================
function x = solveAXEqualsZero(A);
 h  = sym('phi', [1 9]);
 [u,s,v] = svd(A);
% v = v';
 x = v(:,end);
 
%****TO DO **** Write this routine 
