%note the coeffs inside a branch is inside the []

for i=1:length(gp.pop)
    A{i} = tree2evalstr(gp.pop{i}, gp);
end

string2Beval = string(gp.pop);
FD_re = string2Beval;
FD_mu = string2Beval;
coeffies = zeros(length(string2Beval), 1);

string2Beval{11}

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
            
            temp = ""; temp2 = ""; temp_m = "";
            for k = 1:length(str_idx)
                %disp(i)
                coeffies(i,k) = str2num(char_temp(str_idx(k)+1:end_idx(k)-1));
                temp = temp + breakdown(k) + string(coeffies(i,k)) + "+0.000001i";
                temp2 = temp2 + breakdown(k) + string(coeffies(i,k)) + "+0.000001";
                temp_m = temp_m + breakdown(k) + string(coeffies(i,k)) + "+multicomplex([0, 0.000001])";
                
            end
            string2Beval(i) = temp + breakdown(end);
            FD_re(i) = temp2 + breakdown(end);
            FD_mu(i) = temp_m + breakdown(end);
        end
        
        
    %end
end

B = convertStringsToChars(string2Beval);
B2 = convertStringsToChars(FD_mu);
for i=1:length(B)
    C{i} = cellstr(B{i});
    C2{i} = cellstr(B2{i});
end

B2{11}

%% 

gp2 = gp;

gp2.pop = C;

%gp2 = evalfitness(gp2);

h = 0.000001;
Res_cs = imag(gp2.fitness.values)/h; %complex set method
Res_cs(isnan(Res_cs))=0;
Res_cs(isinf(Res_cs))=0;
Res_cs = round(Res_cs, 3);

gp2.pop = C2;
% gp2 = evalfitness(gp2);

Res_fd = (gp2.fitness.values-gp.fitness.values)/h; %check with FD method
Res_fd(isnan(Res_fd))=0;
Res_fd(isinf(Res_fd))=0;
Res_fd = round(Res_fd, 3);
