b= 2;
feature('numCores')
matlabpool('local');

parfor i = 1 : 10000
    o = jsqq(i,b)
end


