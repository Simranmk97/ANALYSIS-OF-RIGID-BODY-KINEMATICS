% Clear all stored variables, clear command screen, set output format to 
% short (i.e., 5 digit format), and close all figures
clearvars
clc
format short 
close all 
 
%% LOAD DATA FILES FOR 3 POSITIONS

pos_1_data = readmatrix('VSKN_Position1.xlsx','Range','B1:D21');
pos_2_data = readmatrix('VSKN_Position2.xlsx','Range','B1:D21');
pos_3_data = readmatrix('VSKN_Position3.xlsx','Range','B1:D21');
% ----------------------------------------------------------------------- % 
 
%% ESTABLISH ORTHONORMAL COORDINATE SYSTEM FOR FEMUR MARKER SET

% Load initial coordinate positions of markers
% position 1
fmt = pos_1_data(1,:); 
fmb = pos_1_data(2,:);   %origin of the markers
fmf = pos_1_data(3,:);
 
% position 2
fmt2 = pos_2_data(1,:); 
fmb2 = pos_2_data(2,:);   
fmf2 = pos_2_data(3,:);
 
% position 3
fmt3 = pos_3_data(1,:); 
fmb3 = pos_3_data(2,:);   
fmf3 = pos_3_data(3,:);

% position 1
a = fmt-fmb; % y axis
b = fmf-fmb; % x axis
vec_y = a;
vec_x = b;
vec_z = cross(vec_x,vec_y);
vec_x = cross(vec_y,vec_z);
e_fmy = (vec_y/norm(vec_y)) %anterior direction
e_fmx = (vec_x/norm(vec_x)) %distal direction
e_fmz = (vec_z/norm(vec_z)) %medial direction
 
% position 2
a2 = fmt2-fmb2; % y axis
b2 = fmf2-fmb2; % x axis
vec_y2 = a2;
vec_x2 = b2;
vec_z2 = cross(vec_x2,vec_y2);
vec_x2 = cross(vec_y2,vec_z2);
e_fmy2 = (vec_y2/norm(vec_y2)) %anterior direction
e_fmx2 = (vec_x2/norm(vec_x2)) %distal direction
e_fmz2 = (vec_z2/norm(vec_z2)) %medial direction
 
% position 3
a3 = fmt3-fmb3; % y axis
b3 = fmf3-fmb3; % x axis
vec_y3 = a3;
vec_x3 = b3;
vec_z3 = cross(vec_x3,vec_y3);
vec_x3 = cross(vec_y3,vec_z3);
e_fmy3 = (vec_y3/norm(vec_y3)) %anterior direction
e_fmx3 = (vec_x3/norm(vec_x3)) %distal direction
e_fmz3 = (vec_z3/norm(vec_z3)) %medial direction
 
% e_fmx,y,z: base vectors derived from femur marker set 
 
% w = cross(vec_x,vec_z)
% w2 = w/norm(w)
% e_fmy
% norm(e_fmy)
% norm(w2) % the unit normal value are the same
 
% ----------------------------------------------------------------------- %
 
%% ESTABLISH ORTHONORMAL COORDINATE SYSTEM FOR TIBIA MARKER SET
% tmt, digitized tibia top - refers to the divot at the top of the L piece
% tmb, digitized tibia bottom - refers to the middle divot on the L piece
% tmf, digitized tibia front - refers to the divot on the L piece closest
%                              to the joint
 
% Load initial coordinate positions of markers
% postion 1
tmt = pos_1_data(4,:); 
tmb = pos_1_data(5,:);   
tmf = pos_1_data(6,:);
 
% position 2
tmt2 = pos_2_data(4,:); 
tmb2 = pos_2_data(5,:);   
tmf2 = pos_2_data(6,:);
 
% position 3
tmt3 = pos_3_data(4,:); 
tmb3 = pos_3_data(5,:);   
tmf3 = pos_3_data(6,:);
 
% Establish an orthonormal coordinate system for the tibia market set.
 
% e_tmx,y,z: base vectors derived from tibia marker set
% postion 1
c = tmt-tmb; % y axis
d = tmf-tmb; % x axis
vec_y2 = c;
vec_x2 = d;
vec_z2 = cross(vec_x2,vec_y2);
vec_x2 = cross(vec_y2,vec_z2);
e_tmy = (vec_y2/norm(vec_y2))
e_tmx = (vec_x2/norm(vec_x2))
e_tmz = vec_z2/norm(vec_z2)
 
% position 2
c2 = tmt2-tmb2; % y axis
d2 = tmf2-tmb2; % x axis
vec_y22 = c2;
vec_x22 = d2;
vec_z22 = cross(vec_x22,vec_y22);
vec_x22 = cross(vec_y22,vec_z22);
e_tmy2 = (vec_y22/norm(vec_y22))
e_tmx2 = (vec_x22/norm(vec_x22))
e_tmz2 = vec_z22/norm(vec_z22)
 
% position 3
c3 = tmt3-tmb3; % y axis
d3 = tmf3-tmb3; % x axis
vec_y23 = c3;
vec_x23 = d3;
vec_z23 = cross(vec_x23,vec_y23);
vec_x23 = cross(vec_y23,vec_z23);
e_tmy3 = (vec_y23/norm(vec_y23))
e_tmx3 = (vec_x23/norm(vec_x23))
e_tmz3 = vec_z23/norm(vec_z23)
% ----------------------------------------------------------------------- % 
 
%% ESTABLISH ORTHONORMAL COORDINATE SYSTEM FOR MOUNTING FRAME CORNERS
% ftc, digitized mounting frame top - divot on upper back corner, near
%                                     femur
% fbc, digitized mounting frame bottom - divot on lower back corner
% ffc, digitized mounting frame front - divot on lower front corner
 
% Load initial coordinate positions of markers
% position 1
ftc = pos_1_data(7,:); 
fbc = pos_1_data(8,:);   
ffc = pos_1_data(9,:);
 
%position 2
ftc2 = pos_2_data(7,:); 
fbc2 = pos_2_data(8,:);   
ffc2 = pos_2_data(9,:);
 
%position 3
ftc3 = pos_3_data(7,:); 
fbc3 = pos_3_data(8,:);   
ffc3 = pos_3_data(9,:);
 

% Establish an orthonormal coordinate system for the mounting framer 
% corners.
 
% e_mfx,y,z: base vectors derived from mounting frame corners
% position 1
e = ftc- fbc; % y axis
f = ffc - fbc; % x axis 
vec_y3 = e;
% vec_x3 = f;
vec_z3 = cross(f,vec_y3);
vec_x3 = cross(vec_y3,vec_z3);
e_mfy = vec_y3/norm(vec_y3)
e_mfx = vec_x3/norm(vec_x3)
e_mfz = vec_z3/norm(vec_z3)
 
%  position 2
 
e2 = ftc2- fbc2; % y axis
f2 = ffc2 - fbc2; % x axis 
vec_y32 = e2;
% vec_x3 = f;
vec_z32 = cross(f2,vec_y32);
vec_x32 = cross(vec_y32,vec_z32);
e_mfy2 = vec_y32/norm(vec_y32)
e_mfx2 = vec_x32/norm(vec_x32)
e_mfz2 = vec_z32/norm(vec_z32)
 
% position 3
 
e3 = ftc3- fbc3; % y axis
f3 = ffc3 - fbc3; % x axis 
vec_y33 = e3;
% vec_x3 = f;
vec_z33 = cross(f3,vec_y33);
vec_x33 = cross(vec_y33,vec_z33);
e_mfy3 = vec_y33/norm(vec_y33)
e_mfx3 = vec_x33/norm(vec_x33)
e_mfz3 = vec_z33/norm(vec_z33)
 
% ----------------------------------------------------------------------- % 
 
%% ESTABLISH EMBEDDED ORTHONORMAL COORDINATE SYSTEM FOR DISTAL FEMUR
% efdm, embedded femur distal medial - divot on medial side of distal femur 
%                                      (side with market sets)
% efdl, embedded femur distal lateral - divot on lateral side of distal
%                                       femur (nearest mounting frame)
% efdp, embedded femor distal posterior - divot on posterior (top) side of 
%                                         distal femur
% efda, embedded femor dital anterior - divot on anterior (bottom) side of 
%                                       distal femur
% efpp, embedded femur proximal posterior - divot on posterior (top) side 
%                                           of proximal femor
% efpa, embedded femur proximal anterior - divot on anterior (bottom) side 
%                                          of proximal femor
% efdml, embedded femur distal medial-lateral mid-point - location of point
%                                         between medial and lateral divots
%                                         within the distal femur
 
% Load initial coordinate positions of markers
 
% position 1
efdm = pos_1_data(10,:); % medial (x direction)
efdl = pos_1_data(11,:); % lateral (x direction)
efdp = pos_1_data(12,:); 
efda = pos_1_data(13,:);
efpp = pos_1_data(14,:);
efpa = pos_1_data(15,:);
 
e_femur_distal = mean([efdp;efda])
 
e_femur_proximal = mean([efpp;efpa])
 
%position 2
 
efdm2 = pos_2_data(10,:); % medial (x direction)
efdl2 = pos_2_data(11,:); % lateral (x direction)
efdp2 = pos_2_data(12,:); 
efda2 = pos_2_data(13,:);
efpp2 = pos_2_data(14,:);
efpa2 = pos_2_data(15,:);
 
e_femur_distal2 = mean([efdp2;efda2])
 
e_femur_proximal2 = mean([efpp2;efpa2])
 
% position 3
efdm3 = pos_3_data(10,:); % medial (x direction)
efdl3 = pos_3_data(11,:); % lateral (x direction)
efdp3 = pos_3_data(12,:); 
efda3 = pos_3_data(13,:);
efpp3 = pos_3_data(14,:);
efpa3 = pos_3_data(15,:);
 
e_femur_distal3 = mean([efdp3;efda3])
 
e_femur_proximal3 = mean([efpp3;efpa3])
 
% Establish an orthonormal coordinate system embedded in the femur. Note 
% that the origin should be inside the part and understand the correct 
% orientation of the axes (see Figure 3 in lab introduction).
 
% e_efx,y,z: embedded base vectors from distal femur
 
% position 1
 
g = e_femur_proximal-e_femur_distal; % z direction
h = efdm-efdl; % x direction
% k = efda-e_femur_distal; % y direction
vec_y4 = cross(g,h);
% vec_z4 = g;
% vec_y4 = k;
vec_x4 = h;
vec_z4 = cross(vec_x4,vec_y4)
e_efx = vec_x4/norm(vec_x4)
e_efy = vec_y4/norm(vec_y4)
e_efz = vec_z4/norm(vec_z4)
 
% position 2
g2 = e_femur_proximal2-e_femur_distal2; % z direction
h2 = efdm2-efdl2; % x direction
% k = efda-e_femur_distal; % y direction
vec_y42 = cross(g2,h2);
% vec_z4 = g;
% vec_y4 = k;
vec_x42 = h2;
vec_z42 = cross(vec_x42,vec_y42)
e_efx2 = vec_x42/norm(vec_x42)
e_efy2 = vec_y42/norm(vec_y42)
e_efz2 = vec_z42/norm(vec_z42)
 
% position 3
g3 = e_femur_proximal3-e_femur_distal3; % z direction
h3 = efdm3-efdl3; % x direction
% k = efda-e_femur_distal; % y direction
vec_y43 = cross(g3,h3);
% vec_z4 = g;
% vec_y4 = k;
vec_x43 = h3;
vec_z43 = cross(vec_x43,vec_y43)
e_efx3 = vec_x43/norm(vec_x43)
e_efy3 = vec_y43/norm(vec_y43)
e_efz3 = vec_z43/norm(vec_z43)
 
 
% ----------------------------------------------------------------------- % 
 
%% ESTABLISH EMBEDDED ORTHONORMAL COORDINATE SYSTEM FOR PROXIMAL TIBIA
% etpm, embedded tibia proximal medial - divot on medial side of proximal
%                                        stibia (side with market sets)
% etpl, embedded tibia proximal lateral - divot on lateral side of 
%                                         proximal tibia (nearest 
%                                         mounting frame)
% etpp, embedded tibia proximal posterior - divot on posterior (top) side 
%                                           of proximal femur
% etpa, embedded tibia proximal anterior - divot on anterior (bottom) side 
%                                          of proximal femur
% etdp, embedded tibia distal posterior - divot on posterior (top) side 
%                                         of distal femur
% etda, embedded tibia distal anterior - divot on anterior (bottom) side 
%                                        of distal femur
 
% Load initial coordinate positions of markers
%position 1
etpm = pos_1_data(16,:); % x direction
etpl = pos_1_data(17,:); % x direction
etpp = pos_1_data(18,:);
etpa = pos_1_data(19,:); % z direction
etdp = pos_1_data(20,:);
etda = pos_1_data(21,:); % z direction
 
e_tibia_proximal = mean([etpp;etpa]);
 
e_tibia_distal = mean([etdp;etda]);
 
%position 2
etpm2 = pos_2_data(16,:); % x direction
etpl2 = pos_2_data(17,:); % x direction
etpp2 = pos_2_data(18,:);
etpa2 = pos_2_data(19,:); % z direction
etdp2 = pos_2_data(20,:);
etda2 = pos_2_data(21,:); % z direction
 
e_tibia_proximal2 = mean([etpp2;etpa2]);
 
e_tibia_distal2 = mean([etdp2;etda2]);
 
%position 3
 
etpm3 = pos_3_data(16,:); % x direction
etpl3 = pos_3_data(17,:); % x direction
etpp3 = pos_3_data(18,:);
etpa3 = pos_3_data(19,:); % z direction
etdp3 = pos_3_data(20,:);
etda3 = pos_3_data(21,:); % z direction
 
e_tibia_proximal3 = mean([etpp3;etpa3]);
 
e_tibia_distal3 = mean([etdp3;etda3]);
 
% Establish an orthonormal coordinate system embedded in the tibia. Note 
% that the origin should be inside the part and understand the correct 
% orientation of the axes (see Figure 4 in lab introduction).
 
% e_etx,y,z: embedded base vectors for proximal tibia
 
i = e_tibia_proximal-e_tibia_distal; % z direction
j = etpm-etpl; % x direction
% l = etpa-e_tibia_proximal; % y direction
% vec_z4 = i;
% vec_x4 = j;
vec_y4 = cross(i,j);
% vec_x4 = cross(vec_y4,i);
vec_z4 = cross(j,vec_y4);
 
e_etx = vec_x4/norm(vec_x4)
e_ety = vec_y4/norm(vec_y4)
e_etz = vec_z4/norm(vec_z4)
 
% position 2
 
i2 = e_tibia_proximal2-e_tibia_distal2; % z direction
j2 = etpm2-etpl2; % x direction
% l = etpa-e_tibia_proximal; % y direction
% vec_z4 = i;
% vec_x4 = j;
vec_y42 = cross(i2,j2);
% vec_x4 = cross(vec_y4,i);
vec_z42 = cross(j2,vec_y42);
 
e_etx2 = vec_x42/norm(vec_x42)
e_ety2 = vec_y42/norm(vec_y42)
e_etz2 = vec_z42/norm(vec_z42)
 
% position 3
 
i3 = e_tibia_proximal3-e_tibia_distal3; % z direction
j3 = etpm3-etpl3; % x direction
% l = etpa-e_tibia_proximal; % y direction
% vec_z4 = i;
% vec_x4 = j;
vec_y43 = cross(i3,j3);
% vec_x4 = cross(vec_y4,i);
vec_z43 = cross(j3,vec_y43);
 
e_etx3 = vec_x43/norm(vec_x43)
e_ety3 = vec_y43/norm(vec_y43)
e_etz3 = vec_z43/norm(vec_z43)
% ----------------------------------------------------------------------- % 
 
%% PLOT MARKER SETS AND DEFINED COORDINATE SYSTEMS
% check the figure to ensure the calculated coordinate systems axes are 
% correct (compare with Figs. 3-5 in lab introduction)
 
% ** check that variable names are correct **
 
 
cs_data = [fmb e_fmx e_fmy e_fmz;
           tmb e_tmx e_tmy e_tmz;
           fbc e_mfx e_mfy e_mfz;
           e_femur_distal e_efx e_efy e_efz;
           e_tibia_proximal e_etx e_ety e_etz];
  
[~] = plot_markers_CS_2021(pos_1_data,cs_data,1);
 
cs_data2 = [fmb2 e_fmx2 e_fmy2 e_fmz2;
           tmb2 e_tmx2 e_tmy2 e_tmz2;
           fbc2 e_mfx2 e_mfy2 e_mfz2;
           e_femur_distal2 e_efx2 e_efy2 e_efz2;
           e_tibia_proximal2 e_etx2 e_ety2 e_etz2];
  
[~] = plot_markers_CS_2021(pos_2_data,cs_data2,2);
 
cs_data3 = [fmb3 e_fmx3 e_fmy3 e_fmz3;
           tmb3 e_tmx3 e_tmy3 e_tmz3;
           fbc3 e_mfx3 e_mfy3 e_mfz3;
           e_femur_distal3 e_efx3 e_efy3 e_efz3;
           e_tibia_proximal3 e_etx3 e_ety3 e_etz3];
  
[~] = plot_markers_CS_2021(pos_3_data,cs_data3,3);
 
%% DETERMINE 4x4 TRANSFORMATION MATRICES
%
%  Use transform_matrix.m function
%     input 1-3: orthonormal axes (e1, e2, e3) that define c.s.
%     transforming from (unprimed c.s.)
%     input 4-6: orthonormal axes (e1, e2, e3) that define c.s.
%     transforming to (primed c.s.)
%     input 7: x,y,z coordinates of c.s. origin transforming from
%     input 8: x,y,z coordiantes of c.s. origin transforming to
 
%% DETERMINE 4x4 TRANSFORMATION MATRIX BETWEEN FEMUR EMBEDDED C.S. and 
%  MARKET SET C.S.
%
%  (Data Analysis Q1)
 
T_fe_fm = transform_matrix(e_efx,e_efy,e_efz,e_fmx,e_fmy,e_fmz,e_femur_distal,fmb)
 
T_fe_fm2 = transform_matrix(e_efx2,e_efy2,e_efz2,e_fmx2,e_fmy2,e_fmz2,e_femur_distal2,fmb2)
 
T_fe_fm3 = transform_matrix(e_efx3,e_efy3,e_efz3,e_fmx3,e_fmy3,e_fmz3,e_femur_distal3,fmb3)
% ----------------------------------------------------------------------- % 
 
%% DETERMINE 4x4 TRANSFORMATION MATRIX BETWEEN THE TIBIA MARKER SET C.S. 
%  AND EMBEDDED C.S.
%
%  (Data Analysis Q2)
 
T_tm_te = transform_matrix(e_tmx, e_tmy, e_tmz, e_etx, e_ety, e_etz,tmb, e_tibia_proximal)
 
T_tm_te2 = transform_matrix(e_tmx2, e_tmy2, e_tmz2, e_etx2, e_ety2, e_etz2,tmb2, e_tibia_proximal2)
 
T_tm_te3 = transform_matrix(e_tmx3, e_tmy3, e_tmz3, e_etx3, e_ety3, e_etz3,tmb3, e_tibia_proximal3)
 
% ----------------------------------------------------------------------- % 
 
%% DETERMINE 4X4 TRANSFORMATION MATRIX BETWEEN FEMUR MARKER SET C.S. AND
%  TIBIA MARKER SET C.S.
%  (Data Analysis Q3)
 
%  ** Must be calculated for all 3 positions of the kinematic knee **
 
% Position 1: ~5 degrees of flexion and 0 degrees of tibial rotation
 
T_fm_tm_1 = transform_matrix(e_fmx,e_fmy,e_fmz,e_tmx,e_tmy,e_tmz,fmb,tmb)
 
% ----------------------------------------------------------------------- % 
% Position 2: ~30 degrees of knee flexion and ~5 degrees tibial rotation
 
T_fm_tm_2 = transform_matrix(e_fmx2,e_fmy2,e_fmz2,e_tmx2,e_tmy2,e_tmz2,fmb2,tmb2)
 
% ----------------------------------------------------------------------- % 

% Position 3: ~90 degrees of knee flexion and ~30 degrees of tibial 
%             rotation
 
T_fm_tm_3 = transform_matrix(e_fmx3,e_fmy3,e_fmz3,e_tmx3,e_tmy3,e_tmz3,fmb3,tmb3)
 
% ----------------------------------------------------------------------- % 
 
%% DETERMINE THE OVERALL 4x4 TRANSFORMATION MATRIX BETWEEN THE FEMUR 
%  EMBEDDED C.S. AND TIBIA EMBEDDED C.S.
%  (Data Analysis Q4)

T_fe_te_1 = T_tm_te*T_fm_tm_1*T_fe_fm
T_fe_te_2 = T_tm_te*T_fm_tm_2*T_fe_fm
T_fe_te_3 = T_tm_te*T_fm_tm_3*T_fe_fm
 
% ----------------------------------------------------------------------- % 
 
%% DETERMINE EULER ANGLES ASSOCIATED WITH MOTION ON THE FEMUR TO TIBIA 
%  EMBEDDED C.S. (i.e., T_fe_te) FOR ALL 3 POSITIONS
%
% naming convention; angles should be calculated in degrees
% psi - CCW rotation about x-axis
% theta - CCW rotation about y-axis
% phi - CCW rotation about z-axis
% see lecture slides when discussing Euler Angles
%
% It is encouraged to calculate the angles in radians and then covert to
% degrees
%

% Position 1:
T1=T_fe_te_1;
theta1=asin(T1(3,1));
psi1=-atan2(T1(3,2)/cos(theta1),T1(3,3)/cos(theta1));
phi1=-atan2(T1(2,1)/cos(theta1),T1(1,1)/cos(theta1));
 
%Conversion
 
theta1 = rad2deg(theta1)
psi1 = rad2deg(psi1)
phi1 = rad2deg(phi1)
 
 
% Position 2:
 
T2=T_fe_te_2;
theta2=asin(T2(3,1));
psi2=-atan2(T2(3,2)/cos(theta2),T2(3,3)/cos(theta2));
phi2=-atan2(T2(2,1)/cos(theta2),T2(1,1)/cos(theta2));
 
%Conversion
 
theta2 = rad2deg(theta2)
psi2 = rad2deg(psi2)
phi2 = rad2deg(phi2)
 
% Position 3:
 
T3=T_fe_te_3;
theta3=asin(T3(3,1));
psi3=-atan2(T3(3,2)/cos(theta3),T3(3,3)/cos(theta3));
phi3=-atan2(T3(2,1)/cos(theta3),T3(1,1)/cos(theta3));
 
%Conversion
 
theta3 = rad2deg(theta3)
psi3 = rad2deg(psi3)
phi3 = rad2deg(phi3)
 
 
% ----------------------------------------------------------------------- %
 
%% CREATE PLOTS FOR CHANGES IN EULER ANGLES ACROSS POSITIONS
% See plotting guidelines in Report Guide document
 
 
figure
position = [1 2 3];
subplot(3,1,1)
psi_a = [psi1, psi2, psi3];
plot(position, psi_a, '-o', 'LineWidth',2)
title('Psi Euler Angles (Plot a)')
xlabel('Position Number')
ylabel('\psi Degrees')
 
subplot(3,1,2)
theta_a = [theta1, theta2, theta3];
plot(position, theta_a, '-o', 'LineWidth',2)
title('Theta Euler Angles (Plot b)')
xlabel('Position Number')
ylabel('\theta Degrees')
 
subplot(3,1,3)
phi_a = [phi1, phi2, phi3];
plot(position, phi_a, '-o', 'LineWidth',2)
title('Phi Euler Angles (Plot c)')
xlabel('Position Number')
ylabel('\phi Degrees')
 
sgtitle('Group 4, Wednesday 2-5pm; Embedded Euler Angle Plots')
 
saveas(gcf,'eulerplots.png')
 
 
%% Calculate the magnitude of the translation vector
 
r = T_fe_te_1(1:3,4);
q = T_fe_te_2(1:3,4);
ll = T_fe_te_3(1:3,4);
 
t1 = norm(r)
t2 = norm(q)
t3 = norm(ll)
 
 
 
%% Manual calculation of angles forr position 1,2,3 (psi and phi)
% psi caluclations for position 1,2,3
 
vector1 = efda-efdp; % femur distal anterior - femur distal posterior 
vector2 = etpa-etpp; % tibia proximal anterior - tibia proximal posterior
 
v1 = norm(vector1);
v2 = norm(vector2);
 
angle1 = rad2deg(acos((dot(vector1,vector2)./(v1*v2))))
 
 
vector3 = efda2-efdp2; % femur distal anterior - femur distal posterior 
vector4 = etpa2-etpp2; % tibia proximal anterior - tibia proximal posterior
 
v3 = norm(vector3);
v4 = norm(vector4);
 
angle2 = rad2deg(acos((dot(vector3,vector4)./(v3*v4))))
 
 
 
vector5 = efda3-efdp3; % femur distal anterior - femur distal posterior 
vector6 = etpa3-etpp3; % tibia proximal anterior - tibia proximal posterior
 
v5 = norm(vector5);
v6 = norm(vector6);
 
angle3 = rad2deg(acos((dot(vector5,vector6)./(v5*v6))))
 
 
% phi caluclations for position 1,2,3
 
vector7 = etda-etdp; % tibia distal anterior - tibia distal posterior 
vector8 = etpa-etpp; % tibia proximal anterior - tibia proximal posterior
 
v7 = norm(vector7);
v8 = norm(vector8);
 
angle4 = rad2deg(acos((dot(vector7,vector8)./(v7*v8))))
 
 
vector9 = etda2-etdp2; % tibia distal anterior - tibia distal posterior 
vector10 = etpa2-etpp2; % tibia proximal anterior - tibia proximal posterior
 
v9 = norm(vector9);
v10 = norm(vector10);
 
angle5 = rad2deg(acos((dot(vector9,vector10)./(v9*v10))))
 
vector11 = etda3-etdp3; % tibia distal anterior - tibia distal posterior 
vector12 = etpa3-etpp3; % tibia proximal anterior - tibia proximal posterior
 
v11 = norm(vector11);
v12 = norm(vector12);
 
angle6 = rad2deg(acos((dot(vector11,vector12)./(v11*v12))))
----------------------------------------------------------------
 
TRANFORMATION_MATRIX.M FUNCTION

function m = transform_matrix(e1,e2,e3,e1p,e2p,e3p,eo,eop)
 
m = [dot(e1p,e1) dot(e1p,e2) dot(e1p,e3) eop(1)-eo(1)
    dot(e2p,e1) dot(e2p,e2) dot(e2p,e3) eop(2)-eo(2)
    dot(e3p,e1) dot(e3p,e2) dot(e3p,e3) eop(3)-eo(3)
    0 0 0 1]
end


 


