function [fitness,gp] = quartic_fitfun(evalstr,gp)
%QUARTIC_FITFUN Fitness function for simple ("naive") symbolic regression on the quartic polynomial y = x + x^2 + x^3 + x^4. 
%   
%   FITNESS = QUARTIC_FITFUN(EVALSTR,GP) returns the FITNESS value of the
%   symbolic expression contained within the cell array EVALSTR using the
%   information in the GP struct.
%   
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
%
%   See also GPDEMO1_CONFIG, GPDEMO1

%extract x and y data from GP struct
% x1 = gp.userdata.x(:,1);
% x2 = gp.userdata.x(:,1);
x =  gp.userdata.xtrain;
pat = 'x(\d+)';
evalstr = regexprep(evalstr,pat,'x(:,$1)'); %$1 makes x1 and x2 things 
y = gp.userdata.ytrain;

%evaluate the tree (assuming only 1 gene is suppled in this case - if the
%user specified multigene config then only the first gene encountered will be used)
try
    eval(['out=' evalstr{1} ';']);
catch
    evalstr{1}
    ans = gp;
end

numData = length(y); % num of data

%fitness is sum of absolute differences between actual and predicted y
%fitness = sum( abs(out-y) );

err = y - out;
fitness = sqrt(((err'*err)/numData));

% fitness = 0;
% for i=1:length(y)
%     fitness = fitness + (out(i) - y(i))^2;
% end
% fitness = sqrt(fitness/numData);

%if this is a post run call to this function then plot graphs
if gp.state.run_completed
    
    figure('name','GPTIPS 2 Multigene regression. Model predictions.','numbertitle','off');
    subplot(1+plotTest+plotValidation,1,1);
    plot(ypredtrain,'Color',[0.85 0.33 0.1]);
    hold on;
    plot(gp.userdata.ytrain,'Color',[0 0.45 0.74]);
    axis tight;
    ylabel('y');
    xlabel('Data point');
    legend('Predicted','Actual');
    title({setname,...
        ['RMS training set error: ' num2str(fitness) '  R^2: ' num2str(r2train)]},'interpreter','tex');
    hold off;
    
end

