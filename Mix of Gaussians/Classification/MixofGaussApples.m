function r =MixGaussApples
close all;
clear all;
tic
 %% Load the Data 

% Note that cells are accessed using curly-brackets {} instead of parentheses ().
Iapples = cell(6,1);
Iapples{1} = 'Apples/Apples_by_kightp_Pat_Knight_flickr.jpg';
Iapples{2} = 'Apples/ApplesAndPears_by_srqpix_ClydeRobinson.jpg';
Iapples{3} = 'Apples/bobbing-for-apples.jpg';
Iapples{6} = 'Apples/Bbr98ad4z0A-ctgXo3gdwu8-original.jpg';
Iapples{5} = 'Apples/Apples_by_MSR_MikeRyan_flickr.jpg';
Iapples{4} = 'Apples/audioworm-QKUJj2wmxuI-original.jpg';
Iapples{7} = 'Apples/apples_red_green.jpg';
Iapples{8} = 'Apples/Apples_and_Carrots.jpg';

IapplesMasks = cell(4,1);
IapplesMasks{1} = 'Masks/Apples_by_kightp_Pat_Knight_flickr.png';
IapplesMasks{2} = 'Masks/ApplesAndPears_by_srqpix_ClydeRobinson.png';
IapplesMasks{3} = 'Masks/bobbing-for-apples.png';
IapplesMasks{4} = 'Masks/Bbr98ad4z0A-ctgXo3gdwu8-original.png';
IapplesMasks{5} = 'Masks/apples_red_green_mask.png';
IapplesMasks{6} = 'Masks/Apples_and_Carrots_Mask.png';


    apple_red = [];
    apple_blue = [];
    apple_green = [];
    %Apple = [apple_red apple_greee apple_blue]';


    Notapple_red = [];
    Notapple_blue = [];
    Notapple_green = [];

for iImage=1:3
    % Could use this index to loop.
    curI = double(imread(  Iapples{iImage}   )) / 255;
    % curI is now a double-precision 3D matrix of size (width x height x 3). 
    % Each of the 3 color channels is now in the range [0.0, 1.0].
    %figure;
    %imagesc(curI)


    curImask = imread(  IapplesMasks{iImage}   );
    % These mask-images are often 3-channel, and contain grayscale values. We
    % would prefer 1-channel and just binary:
    curImask = curImask(:,:,2) > 128;  % Picked green-channel arbitrarily.
    %figure;
    %imshow(curImask)

    % create Red Green and blue pixel vectors.
    red_pixels = curI(:,:,1);
    green_pixels = curI(:,:,2);
    blue_pixels = curI(:,:,3);
    % find non empty pixels and set to apples
    
    apple_red = [apple_red; red_pixels(find(curImask))];
    apple_blue = [apple_blue ;blue_pixels(find(curImask))];
    apple_green = [apple_green; green_pixels(find(curImask))];
    
    % find non empty pixels and set to Notapples
    
    Notapple_red = [Notapple_red; red_pixels(find(not(curImask)))];
    Notapple_blue = [Notapple_blue ;blue_pixels(find(not(curImask)))];
    Notapple_green = [Notapple_green; green_pixels(find(not(curImask)))];
   
end
%Get them in correct form
Apple = [apple_red apple_green apple_blue]';

NotApple = [Notapple_red Notapple_green Notapple_blue]';





%% Learning the Parameters



%fit mixture of Gaussians
%TO DO fill in this routine (below)
%figure;
%  Use four Gaussians for apples
mixGaussEstApple = fitMixGauss(Apple,4);
save('mixGaussEstApple');

lamdaApple = mixGaussEstApple.weight;
meanApple = mixGaussEstApple.mean;
covApple = mixGaussEstApple.cov;

% 6 Gaussians for not apples.
mixGaussEstNotApple = fitMixGauss(NotApple,6);
save('mixGaussEstNotApple');
lamdaNotApple = mixGaussEstNotApple.weight;
meanNotApple = mixGaussEstNotApple.mean;
covNotApple = mixGaussEstNotApple.cov;

%set priors
priorApple = 0.5;
priorNonApple = 0.5;



for I = 1:8
   
     curI = double(imread(  Iapples{I}   )) / 255;
     [imY, imX, imZ] = size( curI);
     
     % In order to use mvnpdf, the matrix needs to get data that it can
     % take, i.e. each row is a data point
     vec_red =  reshape( curI(:,:,1).',[],1);
     vec_green =  reshape( curI(:,:,2).',[],1);
     vec_blue =  reshape( curI(:,:,3).',[],1);
     pixels = [vec_red vec_green vec_blue]';
        %Calulate the liklihood similar to Part A.
      likeapple = lamdaApple(1)*mvnpdf(pixels' ,meanApple(:,1)', covApple(:,:,1) )+lamdaApple(2)*mvnpdf(pixels' ,meanApple(:,2)', covApple(:,:,2) ) +lamdaApple(3)*mvnpdf(pixels' ,meanApple(:,3)', covApple(:,:,3) );

      likenonapple = lamdaNotApple(1)*mvnpdf(pixels' ,meanNotApple(:,1)', covNotApple(:,:,1) )+lamdaNotApple(2)*mvnpdf(pixels' ,meanNotApple(:,2)', covNotApple(:,:,2) ) +lamdaNotApple(3)*mvnpdf(pixels' ,meanNotApple(:,3)', covNotApple(:,:,3) );
     
      % Reshape and store in a cell to be converted to a matrix later on.
      posteriorApple{I} =  vec2mat((likeapple*priorApple)./(likeapple*priorApple + likenonapple*priorNonApple),imX);
     

     
    
end
Image_1 = cell2mat( posteriorApple(1));
Image_2 = cell2mat( posteriorApple(2));
Image_3 = cell2mat( posteriorApple(3));
Image_4 = cell2mat( posteriorApple(4));
Image_5 = cell2mat( posteriorApple(5));
Image_6 = cell2mat( posteriorApple(6));
Image_7 = cell2mat( posteriorApple(7));
Image_8 = cell2mat( posteriorApple(8));
clims = [0, 1];
%% Plotting the posteriori of test.
figure; set(gcf,'Color',[1 1 1]);
subplot(3,2,1); imagesc(double(imread(  Iapples{4}   )) / 255); axis off; axis image;
subplot(3,2,2); imagesc(Image_4, clims); colormap(gray); axis off; axis image;

subplot(3,2,3); imagesc(double(imread(  Iapples{5}   )) / 255); axis off; axis image;
subplot(3,2,4); imagesc(Image_5, clims); colormap(gray); axis off; axis image;

subplot(3,2,5); imagesc(double(imread(  Iapples{6}   )) / 255); axis off; axis image;
subplot(3,2,6); imagesc(Image_6, clims); colormap(gray); axis off; axis image;
saveas( gcf, 'Test_apple_together', 'jpg' )


%% Part D: testing the truth mask.

%Looking at image 4 mask.     
Imask = double(imread(  IapplesMasks{4}   ));
   
 
count = 1;
% loop of a decreasing threshold
Index = length(1:-0.025:0);
for thresh = 1:-0.025:0
    Current_image = Image_6;
    %set all pixels above threshold equal to one.
    Cur_I = Current_image> thresh;
    % TP is where both mask and posterior agree on pixels where apples exist
    TP = (Cur_I == 1 & Imask==1);
    % sum over the number of times this occurs.
    TTP = sum(sum(TP));
    TN = (Cur_I == 0 & Imask==0);
    TTN = sum(sum(TN));
    FP = (Cur_I == 1 & Imask==0);
    TFP = sum(sum(FP));
    FN = (Cur_I == 0 & Imask==1);
    TFN = sum(sum(FN));
    %True Positive rate
    TPR(count) = TTP/(TTP +TFN);
    %Precision
    PPV(count) = TTP/(TTP +TFP);
    %recall
    recall(count) = TPR(count);
    %False Positive rate
    FPR(count) = TFP/(TFP+TFN);
    %F score
    F(count) =  2*PPV(count)*recall(count)/(PPV(count)+recall(count));
    count = count+1;
    
end

% find highest f score and use that to tune your threshold
[F1Score, order]=max(F);

thresh=1 - ((order-1) *0.025);

currentI = Image_6 > thresh;
disp(['Best Threshold: ' num2str(thresh) ', with F score = ' num2str(F1Score)])
fprintf('\n') 
% find TP,TN,FP,FN of the best threshold
    TP = (currentI == 1 & Imask==1);
    TTP = sum(sum(TP));
    disp(['The total number of True Positives is:' num2str(TTP)])
    TN = (currentI == 0 & Imask==0);
    TTN = sum(sum(TN));
    disp(['The total number of True Negatives is:' num2str(TTN)])
    FP = (currentI == 1 & Imask==0);
    TFP = sum(sum(FP));
    disp(['The total number of False Positives is:' num2str(TFP)])
    FN = (currentI == 0 & Imask==1);
    TFN = sum(sum(FN));
    disp(['The total number of False Negatives is:' num2str(TFN)])
  %% D) Plotting images
% output the images.  
figure; set(gcf,'Color',[1 1 1]);
subplot(2,4,1); imagesc(double(imread(  Iapples{6}   )) / 255); axis off; axis image;title('Origninal');
subplot(2,4,2); imshow(IapplesMasks{4}); title('True Mask');
subplot(2,4,3); imagesc(Image_6, clims); colormap(gray);  axis off; axis image;title('Origninal Posterior');
subplot(2,4,4); imagesc(currentI, clims); colormap(gray);  axis off; axis image;title('Best Posterior');
subplot(2,4,5); imshow(TP);title('True Positive');
subplot(2,4,6);imshow(FP);title('False Positive');
subplot(2,4,7);imshow(TN);title('True Negative');
subplot(2,4,8);imshow(FN);title('False Negative');

%Plot the ROC
figure;
plot(FPR,TPR);title('ROC Curve');xlabel('False Positive Rate');ylabel('True Postive Rate');

%Calculate the area under the ROC curve.
AUC = trapz(FPR,TPR);
disp(['The total area under the curve(AUC) is:' num2str(AUC)])   

%% Part F)

% Now do the same for downloaded masks.

 for I_no = 5:6
   Imask = double(imread(  IapplesMasks{I_no}   ))/255;
 
   % since the pixels are slighly less than 1, make them 1.
   Imask = Imask > 0.5;  

    count = 1;
    for thresh = 1:-0.025:0
        Current_image = cell2mat( posteriorApple(I_no+2));
        Cur_I = Current_image> thresh;
    
      TP = (Cur_I == 1 & Imask==1);
      TTP = sum(sum(TP));
      TN = (Cur_I == 0 & Imask==0);
      TTN = sum(sum(TN));
      FP = (Cur_I == 1 & Imask==0);
      TFP = sum(sum(FP));
      FN = (Cur_I == 0 & Imask==1);
      TFN = sum(sum(FN));
    
      TPR(count) = TTP/(TTP +TFN);
      PPV(count) = TTP/(TTP +TFP);
      recall(count) = TPR(count);
      FPR(count) = TFP/(TFP+TFN);
      F(count) =  2*PPV(count)*recall(count)/(PPV(count)+recall(count));
      count = count+1;
    
    end
    [F1Score, order]=max(F);
    thresh=1 - ((order-1) *0.025);

    currentI = cell2mat( posteriorApple(I_no+2)) > thresh;
    disp(['Best Threshold: ' num2str(thresh) ', with F score = ' num2str(F1Score)])
    fprintf('\n') 
    disp(['The total number of True Positives is:' num2str(TTP)])
    TP = (currentI == 1 & Imask==1);
    TTP = sum(sum(TP));
    disp(['The total number of True Positives is:' num2str(TTP)])
    TN = (currentI == 0 & Imask==0);
    TTN = sum(sum(TN));
    disp(['The total number of True Negatives is:' num2str(TTN)])
    FP = (currentI == 1 & Imask==0);
    TFP = sum(sum(FP));
    disp(['The total number of False Positives is:' num2str(TFP)])
    FN = (currentI == 0 & Imask==1);
    TFN = sum(sum(FN));
    disp(['The total number of False Negatives is:' num2str(TFN)])
  %% D) Plotting images
    
    figure; set(gcf,'Color',[1 1 1]);
    subplot(2,4,1); imagesc(double(imread(  Iapples{I_no+2}   )) / 255); axis off; axis image;title('Origninal');
    subplot(2,4,2); imshow(IapplesMasks{I_no}); title('True Mask');
    subplot(2,4,3); imagesc(cell2mat( posteriorApple(I_no+2)), clims); colormap(gray);  axis off; axis image;title('Origninal Posterior');
    subplot(2,4,4); imagesc(currentI, clims); colormap(gray);  axis off; axis image;title('Best Posterior');
    subplot(2,4,5); imshow(TP);title('True Positive');
    subplot(2,4,6);imshow(FP);title('False Positive');
    subplot(2,4,7);imshow(TN);title('True Negative');
    subplot(2,4,8);imshow(FN);title('False Negative');
    

    figure;
    plot(FPR,TPR);title('ROC Curve');xlabel('False Positive Rate');ylabel('True Postive Rate');
    saveas(gcf,'ROC_test3','jpg');
    AUC = trapz(FPR,TPR);
    disp(['The total area under the curve :' num2str(AUC)])
   
 end



function mixGaussEst = fitMixGauss(data,k);
cIter = 1;
[nDim nData] = size(data);

postHidden = zeros(k, nData);

%in the E-M algorithm, we calculate a complete posterior distribution over
%the (nData) hidden variables in the E-Step.  In the M-Step, we
%update the parameters of the Gaussians (mean, cov, w).  

%we will initialize the values to random values
mixGaussEst.d = nDim;
mixGaussEst.k = k;
mixGaussEst.weight = (1/k)*ones(1,k);
% use k-means to initialize faster
[~ ,mixGaussEst.mean] = kmeans(data',k);
mixGaussEst.mean = mixGaussEst.mean';
for (cGauss =1:k)
    mixGaussEst.cov(:,:,cGauss) = (0.5+1.5*rand(1))*eye(nDim,nDim);
end;

                               
updated_loglike = 123456;
int_loglike = 0;
                               
% Changed to while as it made helped to set a tolerence on the Loglikelihood.
% The tolerence is set here to 1.
while abs(updated_loglike - int_loglike) > 1
            int_loglike = updated_loglike;
           %Expectation step
              
                % using the mvnpdf, the for loop for the data points could be removed

                for j = 1:k
                    responsibilities(j,:) = (mixGaussEst.weight(:,j)*mvnpdf(data',...
                    mixGaussEst.mean(:,j)',mixGaussEst.cov(:,:,j)))';
                end
                
                % The bsxfun is a useful operation that enables one to calculate the responsibilities
                % element by element.
                
                % for cData = 1:nData
                %    postHidden(:,cData) = responsibilities(:,cData)./sum(responsibilities(:,cData));
                %end
                % same as above.
                postHidden = bsxfun(@rdivide,responsibilities,sum(responsibilities));
           %Maximization Step

           %loop over all Gaussians
           for cGauss = 1:k
                % Update weighting parameters 
                mixGaussEst.weight(cGauss) = sum(postHidden(cGauss,:))/sum(sum(postHidden)); 
                % Update mean parameters 
                mixGaussEst.mean(:,cGauss) = (postHidden(cGauss,:)*data(:,:)')/sum(postHidden(cGauss,:));
                % Update covarance parameter 
                 mixGaussEst.cov(:,:,cGauss) = zeros(nDim,nDim);
                 
                 % Here the bsxfun is used once again to reduce the number of for loops required.
                 % As you can see below commented out is what was once used, however being able to remove the
                 % resulted in a faster function.
                 p1 = bsxfun(@minus,data,mixGaussEst.mean(:,cGauss));
                 mixGaussEst.cov(:,:,cGauss) = bsxfun(@times,p1, postHidden(cGauss,:))*p1'./sum(postHidden(cGauss,:));
                 % was getting some errors with cov, so had to add a noise term so it would work.
                 mixGaussEst.cov(:,:,cGauss)=((mixGaussEst.cov(:,:,cGauss)+mixGaussEst.cov(:,:,cGauss)')/2)+0.001*eye(nDim);
                 
           end
           %compute Loglikelihood
  
        
        like = zeros(k,nData);
            % here the mvnpdf is made use of again when calculating the log Likelihood.
            % looping over all the Gaussians again.
        for j=1:size(mixGaussEst.weight,2)
                                                
            like(j,:) = mixGaussEst.weight(:,j) * mvnpdf(data',...
                        mixGaussEst.mean(:,j)',mixGaussEst.cov(:,:,j))';
        end
        like = sum(like);        
        updated_loglike = sum(log(like));  
        fprintf('Log Likelihood Iter %d : %4.3f\n',cIter,updated_loglike);
        cIter = cIter +1;
end
   







