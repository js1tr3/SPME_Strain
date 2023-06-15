function up = Up_1(y)
up =  4.3452-1.6518*(y)+1.6225*(y).^2-2.0843*(y).^3+3.5146*(y).^4-2.2166*(y).^5-0.5623*exp(109.451*(y)-100.006);
up = max(up,0); %  to prevent negative Up values when limits are exceeded
end