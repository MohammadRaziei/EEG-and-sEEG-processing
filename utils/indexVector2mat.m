function [x,y,z] = indexVector2mat(index,dim)
if(index > prod(dim)), error('invalid'); end
x = 1 + mod(index-1,dim(1))
y = 1 + mod(floor((index-1)/dim(1)),dim(2))
z = 1 + floor((index-1)/(dim(1)*dim(2)))