% Loic Chappaz                                  Research

function data = preOrbitConvolution(map, gridID, nreturns)

%% Grid definitions

nxdot = map.grids{gridID}.nxdot;
nx    = map.grids{gridID}.nx;  
x1    = map.grids{gridID}.x0;
x1dot = map.grids{gridID}.x0dot;
x2    = map.grids{gridID}.x1;
x2dot = map.grids{gridID}.x1dot;

%% Pre-convolution

data(nxdot,nx) = struct('x',[],'xdot',[],'cross',[]);

l       = 1;
kount   = 1;

for j=1:nx
    for i=1:nxdot        

        ncross      = map.ncross(l); 
        n           = 0;
         
        if ncross
            range       = kount:kount+ncross-1;
            cross_ij    = map.cross(range,:);
 
            % check cross_ij
            x       = cross_ij(:,1);
            xdot    = cross_ij(:,4);
            bc      = x>=x1 & x<=x2 & xdot>=x1dot & xdot<=x2dot;
            cross_ij= cross_ij(bc, [1 4 8]);
            [n, ~]  = size(cross_ij); 
        end          
        
        if n<=nreturns
            cross_ij = [];
        else
            cross_ij = cross_ij(1:nreturns,:);
                
        end
        
        data(i,j).x      = map.X0(l,1);
        data(i,j).xdot   = map.X0(l,4);
        data(i,j).cross  = cross_ij;
        
        l           = l + 1;
        kount       = kount + ncross;
         

        
    end


end
 


