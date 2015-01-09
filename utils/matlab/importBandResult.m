%
% Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
%
function out = importBandResult(fileName)
    fid = fopen(fileName, 'rt');
    out = [];
    while (~feof(fid))
        type = fscanf(fid, '%s[^\n]');
        if strfind(type, 'EK') == 1
            out.EK = scan();
        elseif strfind(type, 'EIGENVECTOR') == 1
            out.DOS = scan();
        end
    end
    
    fclose(fid);
    
    function M = scan()
        Nk = fscanf(fid, '%d[^\n]');
        data = fscanf(fid, '%d %d[^\n]');
        m = data(1);
        n = data(2);
        
        for ik = 1:Nk
            M.k(ik,:) = fscanf(fid, '%f %f %f[^\n]');
            M.M{ik} = zeros(m,n);
            for ii = 1:m
                for jj = 1:n
                    data = fscanf(fid, '%*[ \n\t]%e', 1);
                    M.M{ik}(ii,jj) = data(1);
                end
            end
        end    
    end
end
    
