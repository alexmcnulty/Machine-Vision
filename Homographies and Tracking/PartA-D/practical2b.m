function practical2b
close all;
%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%load in image 
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
             308.9825  236.7646  255.4416  340.7335   281.5895];
         
%define 3D points of plane
XCart = [-50 -50  50  50 0 ;...
          50 -50 -50  50 0;...
           0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',10);
       
%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.

TEst = estimatePlanePose(xImCart,XCart,K);


%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];

%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points

XWireHom = [XWireFrameCart ; ones(1,size(XWireFrameCart,2))];
XWireHom = TEst*XWireHom;
% Normalize
xCamHom = XWireHom(1:end-1, :);
%move points to image coordinates
xImHom = K*xCamHom;
xImCart2 =xImHom(1:2,:)./ repmat(xImHom(3,:),2,1);

figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart2(1,:),xImCart2(2,:),'r.','MarkerSize',10);


figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart2(1,:),xImCart2(2,:),'r.','MarkerSize',10);hold on;

%draw lines between each pair of points

for (cPoint = 1:4)
    %plot a green line between each pair of points
    if cPoint <= 3 
    plot([xImCart2(1,cPoint) xImCart2(1,cPoint+1)],[xImCart2(2,cPoint) xImCart2(2,cPoint+1)],'g-');
    plot([xImCart2(1,cPoint) xImCart2(1,cPoint+4)],[xImCart2(2,cPoint) xImCart2(2,cPoint+4)],'g-');
    plot([xImCart2(1,cPoint+4) xImCart2(1,cPoint+5)],[xImCart2(2,cPoint+4) xImCart2(2,cPoint+5)],'g-');
    else
        plot([xImCart2(1,cPoint) xImCart2(1,cPoint+4)],[xImCart2(2,cPoint) xImCart2(2,cPoint+4)],'g-');
        plot([xImCart2(1,cPoint+1) xImCart2(1,cPoint+4)],[xImCart2(2,cPoint+1) xImCart2(2,cPoint+4)],'g-');
        plot([xImCart2(1,1) xImCart2(1,4)],[xImCart2(2,1) xImCart2(2,cPoint)],'g-');
    end
    %make sure we don't replace with next point
    hold on;
end;
saveas(gcf, 'practical2b', 'jpg')



%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?




function T = estimatePlanePose(xImCart,XCart,K)

%replace this
T = [];

%TO DO Convert Cartesian image points xImCart to homogeneous representation
%xImHom
xImHom = [xImCart ; ones(1,size(xImCart,2))];

%xImHom = K\xImHom;

xCamHom = K\xImHom;
%xCamHom = xCamHom(1:end-1,:);


%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
%xCamHom

%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.
H = calcBestHomography(XCart,xCamHom);

%phiPrime = K\H;

%TO DO Estimate first two columns of rotation matrix R from the first two
%columns of H using the SVD

% SVD in matlab gives you back v not v'!!!!
[u,s,v] = svd(H(:,1:end-1));
s = [1 0; 0 1; 0 0];
R = u*s*v';
 
%TO DO Estimate the third column of the rotation matrix by taking the cross
%product of the first two columns
phi3 = cross(R(:,1), R(:,2));
R = [R phi3];

%TO DO Check that the determinant of the rotation matrix is positive - if
%not then multiply last column by -1.
if det(R)<0
    R(:,3) = -R(:,3)
end
%TO DO Estimate the translation t by finding the appropriate scaling factor k
%and applying it to the third colulmn of H

%phiPrime = K\H;
% D = bsxfun(@rdivide, C, std(A))
lambda = sum(sum(bsxfun(@rdivide,H(:,1:2),R(:,1:2))))/6;
t = H(:,3)/lambda;
if t(3) < 0
    R(:,1) = -R(:,1);
    R(:,2) = -R(:,2);
    t = -t;
end
%TO DO Check whether t_z is negative - if it is then multiply t by -1 and
%the first two columns of R by -1.

%assemble transformation into matrix form
%T  = [R t;0 0 0 1];
T  = [R t;0 0 0 1];
function H = calcBestHomography(pts1Cart, pts2Cart)

%should apply direct linear transform (DLT) algorithm to calculate best
%homography that maps the points in pts1Cart to their corresonding match in 
%pts2Cart

%****TO DO ****: replace this
H = eye(3);
pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))];
pts1Hom(3,:) = [];
% Already converted to homogeneous
pts2Hom = pts2Cart;
%pts2Hom = [pts2Cart; ones(1,size(pts1Cart,2))];
A = [];

for i=1:length(pts2Hom)
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

