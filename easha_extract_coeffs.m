%%
%note the coeffs inside a branch is inside the []

for j=1:length(gp.pop)
    A{j} = tree2evalstr(gp.pop{j}, gp);
end

string2Beval = string(gp.pop); %string to be evaluated
coeffies = zeros(length(string2Beval), 1); % empty array to store the coefficients
FD_re = string2Beval; %string to be evaluated for FD

h = 1e-7; %complex step size
a = 0.1; %step length
fitness = gp.fitness.values(1);
x1 = gp.userdata.x(1);
y = gp.userdata.y(1);

for j=1:length(string2Beval)
    
        char_temp = convertStringsToChars((string2Beval(j))); 
        str_idx = strfind(char_temp, '['); %start index positions of constants in the function
        
        if ~isempty(str_idx)
            end_idx = strfind(char_temp, ']'); %end index positions of constants in the function
            breakdown = ""; % will contain everything except constants
            breakdown(1) = string(char_temp(1:str_idx(1)));
            for k = 1:length(str_idx)-1
                breakdown(k+1) = string(char_temp(end_idx(k):str_idx(k+1)));
            end
            breakdown(length(str_idx)+1) = string(char_temp(end_idx(length(str_idx)):end ));
            
            temp1 = "";
            temp2 = "";
            temp_new = "";
            
            for k = 1:length(str_idx)
                %disp(i)
                coeffies(j,k) = str2num(char_temp(str_idx(k)+1:end_idx(k)-1)); %extracting coefficients 
%                 temp1 = string(char_temp(1:end_idx(k)-1)) + "+1e-7*i" + string(char_temp(end_idx(k):end));
                temp2 = string(char_temp(1:end_idx(k)-1)) + "+multi([0, 1e-2])" + string(char_temp(end_idx(k):end));
                
%                 temp = cellstr(convertStringsToChars(temp1));
%                 evalstr = tree2evalstr(temp,gp);
%                 eval(['out=' evalstr{1} ';']);
%                 
%                 deriv1 = imag(out)/h;
                
                temp = cellstr(convertStringsToChars(temp2));
                evalstr = tree2evalstr(temp,gp);
                eval(['out=' evalstr{1} ';']);
                
                deriv2 = imgn(out)/h;
                
                %p1 = (fitness - deriv1 .* coeffies(j,k))./deriv1;
                p2 = (fitness - deriv2 .* coeffies(j,k))./deriv2;
                %coeff_new1(j,k) = coeffies(j,k) + a*p1;
                coeff_new2(j,k) = coeffies(j,k) + a*p2;
                
                if isinf(coeff_new2(j,k)) || isnan(coeff_new2(j,k))
                    temp_new = temp_new + breakdown(k) + string(coeffies(j,k));
                else
                    temp_new = temp_new + breakdown(k) + string(coeff_new2(j,k));
                end
                
            end
            string2Beval(j) = temp_new + breakdown(end); %contains the new coefficients
            %FD_re(j) = temp2 + breakdown(end);
        end    
    
end

%evaluate fitness with the new coeffs
% keep new constants only for those whose  fitness_new > fitness 

B = convertStringsToChars(string2Beval);
% B2 = convertStringsToChars(FD_re);
for j=1:length(B)
    C{j} = cellstr(B{j});
%     C2{j} = cellstr(B2{j});
end


gp2 = gp;

gp2.pop = C;

gp2.state.run_completed = false;
for i = 1:gp2.runcontrol.pop_size
    
        %preprocess cell array of string expressions into a form that
        %Matlab can evaluate
        evalstr = tree2evalstr(gp2.pop{i},gp);
        
        
        [fitness,gp2] = feval(gp2.fitness.fitfun,evalstr,gp2);
        gp2.fitness.values(i) = fitness;
        
end
gp2.state.run_completed = true;

disp(['          old constant     new constants'])
   for h = 1:length(gp.pop)
      disp( [gp.fitness.values(h), gp2.fitness.values(h)]);
   end
