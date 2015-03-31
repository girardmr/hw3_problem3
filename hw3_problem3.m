clear all;
close all;
clc;

load('anecoicCystData.mat');

data = veraStrct.data;
fs = 20e6;
speed = 1540; %m/s in body
pixel_size_through_depth = 0.5*(speed/fs); 

for ii = 1:max(size(data))
    time_array_all(ii) = ii/fs;
end

for cc = 1:128
for bb = 1:128
    time_array(:,bb,cc) = time_array_all;
end
end

channel = [[-63.5:1:63.5]];

for beam = 1:128
    
for jj = 1:max(size(data)) %jj=row
    
depth = jj*pixel_size_through_depth; %m

data_matrix = data;
[rows_data_matrix col_data_matrix z_data_matrix] = size(data_matrix);

for ii = 1:(length(channel))
    xe(ii) = 0.1953e-3*abs(channel(ii)); 
    d(ii) = (xe(ii)^2+depth^2)^0.5 + depth;
    time_to_point(ii) = d(ii)/speed;
end

delay_matrix(jj,:,beam) = time_to_point; %delays

end

for aa = 1:128
    delayed_channel(1:rows_data_matrix,aa,beam) = interp1(time_array(1:rows_data_matrix,aa,beam),data_matrix(1:rows_data_matrix,aa,beam),delay_matrix(1:rows_data_matrix,aa,beam),'linear');
end


end


for ll = 1:numel(delayed_channel)
    if isnan(delayed_channel(ll))==1
        delayed_channel(ll) = 0;
    end
end

%binary mask
pitch = 0.1953e-3;
f_number = 2
binary_mask = zeros(rows_data_matrix,128);
for r = 1:rows_data_matrix
    focal_length = r*pixel_size_through_depth;
    D = focal_length/f_number;
    chan_width(r) = round(D/pitch);
    chan_width_ind_1(r) = 64-chan_width(r)/2;
    chan_width_ind_2(r) = 64+chan_width(r)/2;
    for c = 1:128
        if c >= chan_width_ind_1(r) && c<= chan_width_ind_2(r)
            binary_mask(r,c) = 1;
        end
    end
end
imagesc(binary_mask);
colormap gray;
title('Binary mask');

binary_mask = repmat(binary_mask,[1,1,128]);
[num_rows num_col num_beams] = size(delayed_channel);

blackman_win = blackman(128);
blackman_win = blackman_win';
apod_blackman = repmat(blackman_win, [num_rows, 1, num_beams]);
data_blackman = delayed_channel.*apod_blackman;
summed_channels_bl = sum(data_blackman,2);
log_compressed_bl = 20*log10(abs(hilbert(summed_channels_bl(:,:))));
figure;
imagesc(log_compressed_bl,[-50 0]);
axis image;
colormap('gray');
title('Without aperture growth');

blackman_win_g = blackman(128);
blackman_win_g = blackman_win_g';
apod_blackman_g = repmat(blackman_win_g, [num_rows, 1, num_beams]);
apod_blackman_g = apod_blackman_g.*binary_mask;
data_blackman_g = delayed_channel.*apod_blackman_g;
summed_channels_bl_g = sum(data_blackman_g,2);
log_compressed_bl_g = 20*log10(abs(hilbert(summed_channels_bl_g(:,:))));
figure;
imagesc(log_compressed_bl_g,[-50 0]);
axis image;
colormap('gray');
title('Aperture growth');



