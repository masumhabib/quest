%
% Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
%
function [X, Y, Z, V] = importPotential(fileName)
    fid = fopen(fileName, 'rt');
    X = [];
    Y = [];
    Z = [];
    V = [];
    
    data = load(fileName);
    X = data(:,1);
    Y = data(:,2);
    Z = data(:,3);
    V = data(:,4);
    
    fclose(fid);    
end