% Paramters
dirname = '/home/kfujii2/Downloads/transform_mp4/';

num_pics = size( dir(strcat(dirname,'*.png')), 1 );
for i=1:1:num_pics
    common_filename = strcat('img',sprintf('%04d',i),'.png');
    src_filename = strcat(dirname,common_filename);
    dst_filename = strcat(dirname, 'dst/', 'dst_', common_filename);
    output_img_size = 40;

    % Load an image
    src_img = imread(src_filename);

    % Convert to a gray image
    tmp_gray = rgb2gray(src_img);

    % Crop the image
    ori_rect_size = min(size(tmp_gray));
    center_x = floor( size(tmp_gray,2)/2 );
    center_y = floor( size(tmp_gray,1)/2 );
    min_x = center_x - floor(ori_rect_size/2);
    min_y = center_y - floor(ori_rect_size/2);
    tmp_croped_img = imcrop(tmp_gray, [min_x, min_y, ori_rect_size, ori_rect_size]);

    % Change resolution
    dst_img = imresize(tmp_croped_img,[output_img_size,output_img_size]);

    % Save the dst_image
    imwrite(dst_img, dst_filename);
end