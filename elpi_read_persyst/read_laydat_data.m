function [hdr, record] = read_laydat_data(lay_file_name, start_sample, samples_to_read)

    hdr = read_laydat_header(lay_file_name);
    
    %-------------------
    %     dat file
    %-------------------
    dat_file_name = strrep(lay_file_name, '.lay', '.dat');
    
    dat_file_ID = fopen(dat_file_name, 'r', 'a');
    
    %read either int32 or short data type
    if (hdr.datatype =='7')
        precision = 'int32';
        precision_num = 4;
    else
        precision = 'int16';
        precision_num = 2;
        %error("Precision unknown")
    end
    
    status = fseek(dat_file_ID, hdr.nChans*precision_num*(start_sample-1), 'bof');
    record_bytes = fread(dat_file_ID,samples_to_read*hdr.nChans, precision);
    fclose(dat_file_ID);
    
    record = reshape(record_bytes, [hdr.nChans, samples_to_read]);
    record = record * hdr.calibration;
end