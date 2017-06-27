function HW2_TrackingAndHomographies


LLs = HW2_Practical9c( 'll' );
LRs = HW2_Practical9c( 'lr' );
ULs = HW2_Practical9c( 'ul' );
URs = HW2_Practical9c( 'ur' );

close all;

% Load frames from the whole video into Imgs{}.
% This is really wasteful of memory, but makes subsequent rendering faster.
LoadVideoFrames

% Coordinates of the known target object (a dark square on a plane) in 3D:
XCart = [-50 -50  50  50;...
          50 -50 -50  50;...
           0   0   0   0];

% These are some approximate intrinsics for this footage.
K = [640  0    320;...
     0    512  256;
     0    0    1];

% Define 3D points of wireframe object.
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];
 
hImg = figure;
       
% ================================================
for iFrame = 1:numFrames
    xImCart = [LLs(iFrame,:)' ULs(iFrame,:)' URs(iFrame,:)' LRs(iFrame,:)'];
    xImCart = circshift( xImCart, 1);

    % To get a frame from footage 
    im = Imgs{iFrame};

    % Draw image and 2d points
    set(0,'CurrentFigure',hImg);
    set(gcf,'Color',[1 1 1]);
    imshow(im); axis off; axis image; hold on;
    plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',15);
    hold off;
    drawnow
    %TO DO Use your routine to calculate TEst the extrinsic matrix relating the
    %plane position to the camera position.
    TEst = estimatePlanePose(xImCart, XCart, K);



    %TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
    %through the projective camera, and drawing lines betweeen the 
    %resulting 2d image points

    XWireHom = [XWireFrameCart ; ones(1,size(XWireFrameCart,2))];
    XWireHom = TEst*XWireHom;
    % Normalize
    xCamHom = XWireHom(1:end-1, :);
    %move points to image coordinates
    xImHom = K*xCamHom;
    xImCart2 =xImHom(1:2,:)./ repmat(xImHom(3,:),2,1);
    
    % TO DO: Draw a wire frame cube using data XWireFrameCart. You need to
    % 1) project the vertices of a 3D cube through the projective camera;
    % 2) draw lines betweeen the resulting 2d image points.
    % Note: CONDUCT YOUR CODE FOR DRAWING XWireFrameCart HERE
    

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
hold off;
drawnow
    
%     Optional code to save out figure
  %  pngFileName = sprintf( '%s_%.5d.png', 'myOutputH', iFrame );
   %  saveas(gcf, pngFileName, 'png');

    
end % End of loop over all frames.
% ================================================

% TO DO: QUESTIONS TO THINK ABOUT...

% Q: Do the results look realistic?
% If not then what factors do you think might be causing this


% Some images don't look well. Need to look at 
% 122, 105, 99 - unfocused, 92
% 83 is crazy, 82-60 not good, 60- 53 top outside image,
% similar with 30's and 40's


% TO DO: your routines for computing a homography and extracting a 
% valid rotation and translation GO HERE. Tips:
%
% - you may define functions for T and H matrices respectively.
% - you may need to turn the points into homogeneous form before any other
% computation. 
% - you may need to solve a linear system in Ah = 0 form. Write your own
% routines or using the MATLAB builtin function 'svd'. 
% - you may apply the direct linear transform (DLT) algorithm to recover the
% best homography H.
% - you may explain what & why you did in the report.

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


