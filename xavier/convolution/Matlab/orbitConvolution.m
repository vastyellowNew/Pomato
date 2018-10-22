% Loic Chappaz                                  Research

function [rgbIM, X, Xdot] = orbitConvolution(data, npass)

%% Grid definitions

[nxdot, nx] = size(data); 

dx    = diff([data(1,1).x data(1,2).x]);
dxdot = diff([data(1,1).xdot data(2,1).xdot]);
x1    = data(1,1).x;
x2    = data(1,end).x;
x1dot = data(1,1).xdot;
x2dot = data(end,1).xdot;
X     = linspace(x1,x2, nx);
Xdot  = linspace(x1dot,x2dot, nxdot);


%% HERE WE GO

for pass = 1:npass

    fprintf('Pass %d\n',pass)
    
    % Initialize with white noise grid if this is the first pass
    if pass==1
        texture0 = random('uniform', 0, 1, 3, nxdot, nx);
% texture0 = random('unid', 255, 3, nxdot, nx);
        rgb      = zeros(size(texture0));
%         rgb      = texture0;
    else
        texture0 = rgb;
         
    end


%% Run convolution


%parpool(4);
%parfor i=1:nics

for i=1:nxdot
    for j=1:nx
        
        if isempty(data(i,j).cross)
            rgb(:,i,j) = [1 1 1]; 
        else
            x        = data(i,j).cross(:,1);
            xdot     = data(i,j).cross(:,2);
            idx      = floor((x-x1)/dx)+1;
            idxdot   = floor((xdot-x1dot)/dxdot)+1;
            nreturns = length(x);
            rgb_ij   = [0 0 0]';
            
            for ii=1:nreturns
                 rgb_ij = rgb_ij + texture0(:,idxdot(ii),idx(ii));   
            end  
 
            rgb(:,i,j) = rgb_ij/nreturns ; 
        end
    end
end

% Plotting
rgbIM = permute(rgb, [2 3 1 ]); 

% figure();  title(['After ', num2str(nreturns),' return'])
% hold on; xlabel('x, nondimensional'); ylabel('x dot, nondimensional')
% image(X, Xdot, rgbIM)

%%  High pass filter

if pass<npass

    disp('Running GIMP unsharp mask')
    imname  = 'im.png';
    imwrite(rgbIM, imname)

    para    = [1 10 0]; 
    mycmnd  = ['gimp -i -b ''(simple-unsharp-mask "',imname,'" ',num2str(para),')'' -b ''(gimp-quit 0)'''];
    system(mycmnd);
    % [imFilt, map]  = imread(imname);imFilt = ind2rgb(imFilt, map);
    rgbIM  = double(imread(imname))/255;


%     figure(); title(['High passed - After ', num2str(nreturns),' return'])
%     hold on; xlabel('x, nondimensional'); ylabel('x dot, nondimensional')
%     image(X, Xdot, rgbIM)

    rgb = permute(rgbIM, [3 1 2]); 

end
 
% pause
end

 


