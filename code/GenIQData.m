function GenIQData(x, filename, fc, fs)

fileID = fopen(filename,'w');
for i = 1:length(x)
    if (i == length(x))
        fprintf(fileID,'%d,%d',real(x(i))*8192,imag(x(i))*8192);
    else
        fprintf(fileID,'%d,%d,',real(x(i))*8192,imag(x(i))*8192);
    end
end
fclose(fileID);

% EOF