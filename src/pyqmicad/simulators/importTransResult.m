%
% Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
%
function out = importTransResult(fileName)
    fid = fopen(fileName, 'rt');
    out = [];
    while (~feof(fid))
        type = fscanf(fid, '%s[^\n]');
        if strfind(type, 'TRANSMISSION') == 1
            out.TE = scan();
        elseif strfind(type, 'CURRENT') == 1
            ib = sscanf(type, 'CURRENT%d');
            out.IE{ib+1} = scan();
        elseif strfind(type, 'DOS') == 1
            out.DOS = scan();
        elseif strfind(type, 'n') == 1
            out.n = scan();
        end
    end
    
    fclose(fid);
    
    function M = scan()
        NE = fscanf(fid, '%d[^\n]');
        N = fscanf(fid, '%d[^\n]');
        for iE = 1:NE
            M.E(iE) = fscanf(fid, '%f[^\n]');
            M.M{iE} = zeros(N,N);
            for ii = 1:N
                for jj = 1:N
                    data = fscanf(fid, '%*[ \n\t](%e,%e)', 2);
                    M.M{iE}(ii,jj) = data(1) + 1i*data(2);
                end
            end
        end    
    end
end
    
