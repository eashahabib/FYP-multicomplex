%clc; clear; 
close all;

%gp = rungp(@spatial_evol_config);
gp = rungp(@gpdemo1_config)

popbrowser(gp, 'train')

runtree(gp, 'best')


gppretty(gp, 'best')

% modelStruct = gpmodel2struct(gp, 'best')
% 
% modelStruct.train
% 
% %visualing
% modelSym = gpmodel2sym(gp, 'best')
% vpa(modelSym, 2)
% 
% h = ezsurf(modelSym, -5:5, 300);
% h.LineStyle = 'none';
% colormap winter;
% light;
% material shiny;
% h.FaceAlpha = 0.5;
% 
% hold; plot3(gp.userdata.xtest(:,1), gp.userdata.xtest(:,2), gp.userdata.ytest, 'mx');
% 
% figure;
% spatialEvoFunc = sym('1/(1+x1^-4) + 1/(1+x2^-4)')
% modelDiff = spatialEvoFunc - modelSym;
% h2 = ezsurf(modelDiff, -5:5, 300);
% h2.LineStyle = 'none';
% colormap winter;
% light;
% material shiny;
