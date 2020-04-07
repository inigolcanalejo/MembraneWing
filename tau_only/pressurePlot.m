clc
clear

d=importdata('xp');
x=d(:,2);
p=d(:,3);

scatter(x,p)