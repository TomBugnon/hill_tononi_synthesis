clear
close all

% Parameters
dirname = '../movie/';
dst_png_folder = 'dst/';
scrbl_xy_folder = 'scrbl_xy/';
scrbl_t_folder = 'scrbl_t/';
prefix = 'moving_rectangle2';
num_pixel = 40;
output_frame_num = 100;
num_spikes = 200;
target = 'white';

frame_idx = 0;
rand_t_idx = randperm(output_frame_num);
counts = zeros(output_frame_num,1);


% Create folder
root_folder_name = strcat(dirname,'/',prefix, '_rate', num2str(num_spikes));
dst_folder_name = strcat(dirname,'/',prefix, '_rate', num2str(num_spikes),'/dst/');
scrbl_xy_folder_name = strcat(dirname,'/',prefix, '_rate', num2str(num_spikes),'/scrbl_xy/');
scrbl_t_folder_name = strcat(dirname,'/',prefix, '_rate', num2str(num_spikes),'/scrbl_t/');


if ~exist(root_folder_name)
    mkdir(root_folder_name)
end
if ~exist(dst_folder_name)
    mkdir(dst_folder_name)
end
if ~exist(scrbl_xy_folder_name)
    mkdir(scrbl_xy_folder_name)
end
if ~exist(scrbl_t_folder_name)
    mkdir(scrbl_t_folder_name)
end


for i=1:1:output_frame_num
    frame_idx = frame_idx + 1; 
    try 
        % Load image
        src_filename = strcat(dirname,prefix,'/',prefix, '_', sprintf('%02d',frame_idx),'.png');
        I = imread(src_filename);
        % Display image
        %imshow(I)
        title( strcat('frame_idx=',num2str(frame_idx), ' i=',num2str(i)) )
    catch
        frame_idx = 1;
    end

    %{
    % Crop image to be a rectangle images
    src_rect_size = min( size(I) );
    x_center = floor( size(I,2) / 2 );
    y_center = floor( size(I,1) / 2 );
    x_min = x_center - floor(src_rect_size/2);
    y_min = y_center - floor(src_rect_size/2);
    tmp_croped_img = imcrop(I,map, [x_min, y_min, src_rect_size, src_rect_size]);

    %imshow(tmp_croped_img, map);
    %title( strcat('frame_idx=',num2str(frame_idx), ' i=',num2str(i)) )
    %pause(0.1)

    % Convert image size
    small_img = imresize(tmp_croped_img, [num_pixel, num_pixel]);
    %}
    small_img = rgb2gray(I);
    %imshow(small_img)
    
    % Convert intensity data
    % --- "map" represents [R G B], so if map = [1 1 1], that is black.
    % --- So, after the following precessing,
    % --- high value pixels(max=255) in "dst_img" correspond black pixels in "small_img". 
    dst_img = zeros(size(small_img));
    dst_img = small_img;
    %{
    intensity_map = floor( 255*mean(map,2) );
    for j=1:1:length(intensity_map)
        idx = find(small_img==j);
        dst_img(idx) = intensity_map(j);
    end
    %}
    % Limit the number of spikes by choosing pixels randomly
    fired_threshold = 255*(4/5); % decided by no specific reason
    %fired_threshold = 250; % decided by no specific reason
    fired_pixel = find(dst_img > fired_threshold);
    if length(fired_pixel)>num_spikes
        tmp = randperm( length(fired_pixel), (length(fired_pixel)-num_spikes) );
        ignored_pixel = fired_pixel(tmp);
        dst_img(ignored_pixel) = 0; % turn off
    end
        
    imshow(dst_img);
    title( strcat('frame_idx=',num2str(frame_idx), ' i=',num2str(i)) )
    pause(0.1)
    
    
    %--- Save png image
    dst_filename = strcat(dst_folder_name, prefix, '_',sprintf('%04d',i), '.png');
    imwrite(dst_img, dst_filename)
    
    %--- Scramble spatial information
    % scrbl_xy_img = zeros(size(dst_img));
    rand_xy_idx = randperm(numel(dst_img));
    scrbl_xy_img = dst_img(rand_xy_idx);
    scrbl_xy_img = reshape(scrbl_xy_img,[num_pixel,num_pixel]);

    scrbl_xy_filename = strcat(scrbl_xy_folder_name, prefix,'_', sprintf('%04d',i), '.png');
    imwrite(scrbl_xy_img,scrbl_xy_filename);
                           
    %--- Scramble temporal information
    scramble_timestamp = rand_t_idx(i);
    scrbl_t_filename = strcat(scrbl_t_folder_name, prefix,'_', sprintf('%04d',scramble_timestamp), '.png');
    imwrite(dst_img,scrbl_t_filename);
    
    counts(i) = length(find(dst_img>0));
    
end