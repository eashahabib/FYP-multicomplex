function gp = gp_userfcn(gp)
%GP_USERFCN Calls a user defined function once per generation if one has been 
%specified in the field GP.USERDATA.USER_FCN.
% 
%   Remarks:
%
%   The user defined function must accept (as 1st) and return (as 2nd) the
%   GP structure as arguments.
%
%   Example:
%
%   In multigene symbolic regression the function
%   regressmulti_fitfun_validate can be called once per generation to
%   evaluate the best individual so far on a validation ('holdout') data
%   set.
%
%   Copyright (c) 2009-2015 Dominic Searson 
%
%   GPTIPS 2
% 
% C = gp.fitness.returnvalues{gp.state.current_individual}; %current coeffs
% evalstr = tree2evalstr(gp.pop{250},gp);
% 
% y = gp.userdata.ytrain;
% numData = gp.userdata.numytrain; % num of data
% numGenes = numel(evalstr);
% 
% pat = 'x(\d+)';
% evalstr = regexprep(evalstr,pat,'gp.userdata.xtrain(:,$1)+0.000001i')
% 
% 
% %eval each gene in the current individual
% for i = 1:numGenes
%     ind = i + 1;
%     eval(['geneOutputs(:,ind)=' evalstr{i} ';']);
%     geneOutputs(:,ind) = imag(geneOutputs(:,ind))/0.000001;
%     
%     %check for nonsensical answers and break out early with an 'inf' if so
%     if  any(~isfinite(geneOutputs(:,ind))) || any(~isreal(geneOutputs(:,ind)))
%         fitness = Inf;
%         gp.fitness.returnvalues{gp.state.current_individual} = [];
%         return
%     end
% end
% 
% geneOutputs
% 
% %only calc. weighting coeffs during an actual run or if forced
% 
%     %set gp.userdata.bootSample to true to resample data for weights computation
%     
%     %prepare LS matrix
%     
%         goptrans = geneOutputs';
%         prj = goptrans * geneOutputs;
%     
%     %calculate tree weight coeffs using SVD based least squares
%     %normal equation
%     
%             theta = pinv(prj) * goptrans * y;
%         
%         %theta = [];
%         %fitness = Inf;
%         %gp.fitness.returnvalues{gp.state.current_individual} = [];
%         %return;
%     
%   

string2Beval = string(gp.pop);
coeffies = zeros(length(string2Beval), 1);
FD_re = string2Beval;

for i=1:length(string2Beval)
    %for j=1:length(string2Beval{i})
        char_temp = convertStringsToChars((string2Beval(i)));
        str_idx = strfind(char_temp, '[');
        
        if ~isempty(str_idx)
            end_idx = strfind(char_temp, ']');
            breakdown = "";
            breakdown(1) = string(char_temp(1:str_idx(1)));
            for k = 1:length(str_idx)-1
                breakdown(k+1) = string(char_temp(end_idx(k):str_idx(k+1)));
            end
            breakdown(length(str_idx)+1) = string(char_temp(end_idx(length(str_idx)):end ));
            
            temp = ""; 
            temp2 = "";
            for k = 1:length(str_idx)
                %disp(i)
                coeffies(i,k) = str2num(char_temp(str_idx(k)+1:end_idx(k)-1));
                temp = temp + breakdown(k) + string(coeffies(i,k)) + "+0.000001i";
                temp2 = temp2 + breakdown(k) + string(coeffies(i,k)) + "+0.000001";
                
            end
            string2Beval(i) = temp + breakdown(end);
            FD_re(i) = temp2 + breakdown(end);
        end
        
    %end
end

B = convertStringsToChars(string2Beval);
B2 = convertStringsToChars(FD_re);
for i=1:length(B)
    C{i} = cellstr(B{i});
    C2{i} = cellstr(B2{i});
end


gp2 = gp;

gp2.pop = C;

gp2 = evalfitness(gp2);


h = 0.000001;
Res_cs = imag(gp2.fitness.values)/h; %complex set method
Res_cs(isnan(Res_cs))=0;
Res_cs(isinf(Res_cs))=0;
Res_cs = round(Res_cs, 3);



gp2.pop = C2;
gp2 = evalfitness(gp2);

Res_fd = (gp2.fitness.values-gp.fitness.values)/h; %check with FD method
Res_fd(isnan(Res_fd))=0;
Res_fd(isinf(Res_fd))=0;
Res_fd = round(Res_fd, 3);


%NaN happened in complex set when for real values, it was Inf
    
if ~isempty(gp.userdata.user_fcn)
    [~,gp] = feval(gp.userdata.user_fcn,gp);
end

%diff = theta-C