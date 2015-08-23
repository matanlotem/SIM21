function [ zmax ] = getzmax( zc,n )

 zmax=((1+zc)*((1-(n+1)^-2)/(1-n^-2)))-1;
end

