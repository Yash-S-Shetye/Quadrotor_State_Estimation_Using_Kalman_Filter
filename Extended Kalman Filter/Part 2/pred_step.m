function [covarEst,uEst] = pred_step(uPrev,covarPrev,angVel,acc,dt)
%covarPrev and uPrev are the previous mean and covariance respectively
%angVel is the angular velocity
%acc is the acceleration
%dt is the sampling time

pos_x = double(uPrev(1,1));
pos_y = double(uPrev(2,1));
pos_z = double(uPrev(3,1));
qx = double(uPrev(4,1));
qy = double(uPrev(5,1));
qz = double(uPrev(6,1));
pdot_x = double(uPrev(7,1));
pdot_y = double(uPrev(8,1));
pdot_z = double(uPrev(9,1));
wm_x = double(angVel(1,1));
wm_y = double(angVel(2,1));
wm_z = double(angVel(3,1));
am_x = double(acc(1,1));
am_y = double(acc(2,1));
am_z = double(acc(3,1));

R = [cos(qy)*cos(qz), cos(qz)*sin(qx)*sin(qy) - cos(qx)*sin(qz), sin(qx)*sin(qz) + cos(qx)*cos(qz)*sin(qy);
    cos(qy)*sin(qz), cos(qx)*cos(qz) + sin(qx)*sin(qy)*sin(qz), cos(qx)*sin(qy)*sin(qz) - cos(qz)*sin(qx);
    -sin(qy), cos(qy)*sin(qx), cos(qx)*cos(qy)];
G = pinv([cos(qy)*cos(qz), -sin(qz), 0;
    cos(qy)*sin(qz), cos(qz), 0;
    -sin(qy), 0, 1])* R;


%%

ngx = 0.001;
ngy = 0.001;
ngz = 0.001;
nax = 0.001;
nay = 0.001;
naz = 0.001;
bgx = 0.001;
bgy = 0.001;
bgz = 0.001;
bax = 0.001;
bay = 0.001;
baz = 0.001; %found by tuning

ng = [ngx; ngy; ngz];
na = [nax; nay; naz];
nbg = [bgx; bgy; bgz];
nba = [bax; bay; baz];
%%
X3 = [uPrev(7,1); uPrev(8,1); uPrev(9,1)];
X4 = [uPrev(10,1); uPrev(11,1); uPrev(12,1)];
X5 = [uPrev(13,1); uPrev(14,1); uPrev(15,1)];  
wm = [angVel(1,1); angVel(2,1); angVel(3,1)];
am = [acc(1,1); acc(2,1); acc(3,1)];

Xdot = [X3; (G*(wm - X4 -ng)); [0; 0; -9.81]+(R*(am-X5-na)); nbg; nba];

%%
At =[0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, cos(qz)*sin(qy)*(bgx + ngx - wm_x), cos(qz)*(bgy + ngy - wm_y) + cos(qy)*sin(qz)*(bgx + ngx - wm_x), 0, 0, 0, -cos(qy)*cos(qz), sin(qz), 0, 0, 0, 0;
    0, 0, 0, 0, sin(qy)*sin(qz)*(bgx + ngx - wm_x), sin(qz)*(bgy + ngy - wm_y) - cos(qy)*cos(qz)*(bgx + ngx - wm_x), 0, 0, 0, -cos(qy)*sin(qz), -cos(qz), 0, 0, 0, 0;
    0, 0, 0, 0, cos(qy)*(bgx + ngx - wm_x), 0, 0, 0, 0, sin(qy), 0, -1, 0, 0, 0;
    0, 0, 0, - (sin(qx)*sin(qz) + cos(qx)*cos(qz)*sin(qy))*(bay - am_y + nay) - (cos(qx)*sin(qz) - cos(qz)*sin(qx)*sin(qy))*(baz - am_z + naz), cos(qz)*sin(qy)*(bax - am_x + nax) - cos(qx)*cos(qy)*cos(qz)*(baz - am_z + naz) - cos(qy)*cos(qz)*sin(qx)*(bay - am_y + nay), (cos(qx)*cos(qz) + sin(qx)*sin(qy)*sin(qz))*(bay - am_y + nay) - (cos(qz)*sin(qx) - cos(qx)*sin(qy)*sin(qz))*(baz - am_z + naz) + cos(qy)*sin(qz)*(bax - am_x + nax), 0, 0, 0, 0, 0, 0, -cos(qy)*cos(qz), (cos(qx)*sin(qz) - cos(qz)*sin(qx)*sin(qy)), - sin(qx)*sin(qz) - cos(qx)*cos(qz)*sin(qy);
    0, 0, 0, (cos(qz)*sin(qx) - cos(qx)*sin(qy)*sin(qz))*(bay - am_y + nay) + (cos(qx)*cos(qz) + sin(qx)*sin(qy)*sin(qz))*(baz - am_z + naz), sin(qy)*sin(qz)*(bax - am_x + nax) - cos(qx)*cos(qy)*sin(qz)*(baz - am_z + naz) - cos(qy)*sin(qx)*sin(qz)*(bay - am_y + nay), (cos(qx)*sin(qz) - cos(qz)*sin(qx)*sin(qy))*(bay - am_y + nay) - (sin(qx)*sin(qz) + cos(qx)*cos(qz)*sin(qy))*(baz - am_z + naz) - cos(qy)*cos(qz)*(bax - am_x + nax), 0, 0, 0, 0, 0, 0, -cos(qy)*sin(qz), - cos(qx)*cos(qz) - sin(qx)*sin(qy)*sin(qz), (cos(qz)*sin(qx) - cos(qx)*sin(qy)*sin(qz));
    0, 0, 0, cos(qy)*sin(qx)*(baz - am_z + naz) - cos(qx)*cos(qy)*(bay - am_y + nay), cos(qy)*(bax - am_x + nax) + cos(qx)*sin(qy)*(baz - am_z + naz) + sin(qx)*sin(qy)*(bay - am_y + nay), 0, 0, 0, 0, 0, 0, 0, sin(qy), -cos(qy)*sin(qx), -cos(qx)*cos(qy);
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];


Ft = eye(15) + dt*At; %Ft for prediction step

%%

Ut = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    -cos(qy)*cos(qz), sin(qz), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    -cos(qy)*sin(qz), -cos(qz), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    sin(qy), 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -cos(qy)*cos(qz), cos(qx)*sin(qz) - cos(qz)*sin(qx)*sin(qy), - sin(qx)*sin(qz) - cos(qx)*cos(qz)*sin(qy), 0, 0, 0, 0, 0, 0;
    0, 0, 0, -cos(qy)*sin(qz), - cos(qx)*cos(qz) - sin(qx)*sin(qy)*sin(qz), cos(qz)*sin(qx) - cos(qx)*sin(qy)*sin(qz), 0, 0, 0, 0, 0, 0;
    0, 0, 0, sin(qy), -cos(qy)*sin(qx), -cos(qx)*cos(qy), 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1];

Qd = diag([ngx ngy ngz nax nay naz bgx bgy bgz bax bay baz])*dt; %Qd matrix found by tuning the values

%%

uEst = double(uPrev + Xdot*dt); %Estimated mean
covarEst = double(Ft*covarPrev*transpose(Ft)+ Ut*Qd*transpose(Ut)); %Estimated Covariance
    
end

