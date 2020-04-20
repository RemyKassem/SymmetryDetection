clear all; rng(10000); %clear memory and set seed
%% Asymmetric Images
size = 270;
I=zeros(size,size);
I(60:180,90:120)=.3;
I(150:180,120:150)=.3;
I(120:180,150:180)=.3;
K = imgaussfilt(I,10);
figure(1), imshow(I)
title('Asymmetric shape')

%% Create C2 images Part 1
% (If size = 400) Define I such that x is in (50,200) and y is in (50,350) 

% Data 1 ) 
% size = 400;
% I = zeros(size,size); 
% I(  50:80,60:90) = .5;
% I(100:200,115:350) = .3;

% Data 2 
% size = 400;
% I = zeros(size,size); 
% I(60:130,70:270)=.5;
% I(150:190,300:330)=.2;
% I(80:170,70:100) = 1;

% Data 3  (4 F's with two bright spots) 
% size = 270;
% I = zeros(size,size);
% I(40:120,145:175)=.3; 
% I(40:60,175:205)=.7;
% I(80:100,175:205)=.3; 
% I(40:120,95:125)=.3; 
% I(40:60,65:95)=.3;
% I(80:100,65:95)=.3; 


%% Create C2 images Part 2
% copy over to other side

% J = imrotate(I,180);
% K = I +J; 
% K = imgaussfilt(K,10); % gaussian smoothing
% figure(1), imshow(I)
% title('Starting image')
% figure(2), imshow(K)
% title('C2 symmetry')

%% Create D2 images Part1
% (If size = 400) Define I such that x is in (50,200) and y is in (50,200)

% Data 1 
% size = 400;
% I = zeros(size,size); 
% I(50:80,60:90) = .5;
% I(70:180,115:150) = .3;

% Data 2 
% size = 400; 
% I = zeros(size,size); 
% I(60:90,90:200)=.7;
% I(120:180,60:140)=.3;
% I(100:120,120:140)=1;


% Data 3 (4 Fs) 
% size = 270; 
% I = zeros(size,size); 
% I(40:120,95:125)=.3; 
% I(40:60,65:95)=.3;
% I(80:100,65:95)=.3; 

% Data 4 (4 Fs with 4 bright spots) 
% size = 270;
% I = zeros(size,size); 
% I(40:120,95:125)=.3; 
% I(40:60,65:95)=.7;
% I(80:100,65:95)=.3; 


%% Create D2 images Part 2
% copy over to other quadrants
% J = I +flip(I,2);
% K = J + imrotate(J, 180);
%K = imgaussfilt(K,10); % gaussian smoothing
% figure(1), imshow(I)
% title('Starting image')
% figure(2), imshow(K)
%title('D2 symmetry')


%% Create C4 Images Part I
% Data 1
% size = 400;
% I = zeros(size,size); 
% I(50:80,60:90) = .7;
% I(70:180,115:150) = .3;
%% Create C4 Images Part II

% K = I + imrotate(I, 90) +imrotate(I, 180) + imrotate(I, 270) ;
% K = imgaussfilt(K,10); % gaussian smoothing
% figure(1), imshow(I)
% title('I') 
% figure(2), imshow(K)
%title('C4 symmetry') 

%% Create D1 Images Part I
%Define I such that x is in (50,350) and y is in (50,200)
% Data 1 
% size = 400;
% I = zeros(size,size); 
% I(50:350,60:80) = .2;
% I(70:120,150:190) = .7;
% I(70:110,100:140) = .5;
%% Create D1 Images Part II
% K = I + flip(I,2);
% K = imgaussfilt(K,10);
% figure(1), imshow(I)
% title('I') 
% figure(2), imshow(K)
% title('D1 symmetry') 

%% Sinograms
n=5000;
theta= rand(n,1)*360;

theta = sort(theta,'ascend'); 

M = zeros(n,size); 
for k = 1:n %each row of M is a sinogram
    M(k,:) = sum( imrotate(K,theta(k),'bilinear','crop') );
end

%% construct affintiy matrix (Set-up)
N = 15; %bandwidth parameter
[idx, disX] = knnsearch(M, M, 'k', N);


%% construct affinty matrix W (Self-tuning)
tic
sigma = disX(:,N);
dis = squareform(pdist(M));
toc

tic
W = zeros(n,n);
for i = 1:n
    for j = 1:n
        W(i,j) = exp(-dis(i,j)^2/(2*sigma(i)*sigma(j)));
    end
end
toc
figure(3), imagesc(W)
title('Image of Affinity Matrix')

%% eig( D^{-1}*W )
tic
dW = sum(W,2);
[v,d] = eig(W, diag(dW));
[Lam, ind] = sort(diag(d),'descend');
Psi = v(:,ind); % order the eigenvectors
toc

%% Diffusion Map
%compute top 3 nontrivial eigenvectors
vec1 = Psi(:,2);
vec2 = Psi(:,3);
vec3 = Psi(:,4);

t = 0; %define time t
x = (Lam.^t).*vec1;
y = (Lam.^t).*vec2;
z = (Lam.^t).*vec3;


%% Create colored figures for Diffusion Map
c = linspace(1,10,length(x));
figure(4), scatter(x,y,[],c,'o')
title('Diffusion map for Shape 2d')
xlabel('psi_1')
ylabel('psi_2')
figure(5), scatter3(x,y,z,[],c,'o')
title('Diffusion map for Shape 3d')
xlabel('psi_1')
ylabel('psi_2')
zlabel('psi_3')


