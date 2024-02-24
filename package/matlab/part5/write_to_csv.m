function write_to_csv(C, odor_names, filename)

fid = fopen(filename, 'w+');

for (i=1:length(C))
    animal_name = strcat('F_',num2str(i), ',');   
    for(j=1:size(C{i},1))
        odor_name = strcat(odor_names{j}, ',');  
        fprintf(fid, animal_name);
        fprintf(fid, odor_name);
        for(k=1:size(C{i},2))
            fprintf(fid, '%f', C{i}(j,k));
             if(k<size(C{i},2)) fprintf(fid, ', '); end
        end
       fprintf(fid, '\n');
    end
end

fclose(fid);
