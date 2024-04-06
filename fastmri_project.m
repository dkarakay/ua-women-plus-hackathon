% Deniz Karakay 
% ECE 529 Project

clearvars
close all;
clear all;

% Save the results to a .mat file or not
is_save = true;

% File name for saving the results
file_name = 'test3';

% Create a new folder under 'results/phantom' with the name 'file_name'
folder_name = ['results/fastmri/', file_name];

if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

% Read the original image
filename = 'file1000200.h5';

% Get the filename without the extension
[~, filename_w, ~] = fileparts(filename);

fileinfo = h5info(filename);

% Dataset K-Space
ds = '/kspace';

% Get Slice
slice = 34;
undersampling_rates = [2];
result_index = 1;

results = struct('filename', {}, 'slice', {}, 'undersampling_rate', {}, 'sense_image', {}, 'mse', {}, 'psnr', {}, 'snr', {});

for undersampling_rate = undersampling_rates

    % Get the k-space data of all slices
    data = h5read(filename, ds);

    % Read data as complex numbers
    data_complex = complex(data.r, data.i);

    % Permute the data to match the width, height, coils, slices order
    k_space_data_all = permute(data_complex, [2, 1, 3, 4]);

    % Get the size of the k-space data
    [size_x, size_y, n_coils, n_slices] = size(k_space_data_all);

    % Get the k-space data of given slice
    k_space = k_space_data_all(:, :, :, slice);

    font_title = 16;

    % Plot the k-space data
    f = figure;
    % Set the figure window's size ([left bottom width height])
    f.Position = [0 0 1000, 600];

    % Only 1-6-11 coils
    colormap('jet');
    subplot(1, 3, 1);
    imagesc(abs(log(k_space(:, :, 1))));
    title('K-Space Image: Coil 1', 'FontSize', font_title);
    colorbar();

    subplot(1, 3, 2);
    imagesc(abs(log(k_space(:, :, 6))));
    title('K-Space Image: Coil 6', 'FontSize', font_title);
    colorbar();

    subplot(1, 3, 3);
    imagesc(abs(log(k_space(:, :, 11))));
    title('K-Space Image: Coil 11', 'FontSize', font_title);
    colorbar();

    sgtitle('K-Space Data in Log Scale', 'FontSize', 20);

    % Combined K-Space image plot
    f = figure;
    % Set the figure window's size ([left bottom width height])
    f.Position = [0 0 400, 700];

    colormap('jet');
    imagesc(abs(log(custom.rss(k_space))));
    title('Combined K-Space Image in Log Scale', 'FontSize', font_title);
    colorbar();

    % Get the coil images
    coil_images = custom.ifft2op(k_space);
    % Combine the coil images using the root sum of squares (RSS) method
    img_combined = custom.rss(coil_images);

    % Estimate the coil sensitivity maps
    coil_maps = custom.coil_map_estimation(coil_images, img_combined);

    % Show all coils total 15
    f = figure;
    % Set the figure window's size ([left bottom width height])
    f.Position = [0 0 1000, 600];

    for i = 1:15
        colormap('gray');
        subplot(3, 5, i);
        imagesc(abs(coil_images(:, :, i)));
        title(['Coil Image: ', num2str(i)]);
        colorbar();
    end

    % Show the combined image
    f = figure;
    % Set the figure window's size ([left bottom width height])
    f.Position = [0 0 400, 700];

    colormap('gray');
    imagesc(abs(img_combined));
    title('Combined Coil Images', 'FontSize', font_title);
    colorbar();

    % Show all coil sensitivity maps
    f = figure;
    % Set the figure window's size ([left bottom width height])
    f.Position = [0 0 1000, 600];

    for i = 1:15
        subplot(3, 5, i);
        imagesc(abs(coil_maps(:, :, i)));
        title(['Coil Map: ', num2str(i)]);
        colorbar();
    end

    % Show the combined coil sensitivity map
    f = figure;
    % Set the figure window's size ([left bottom width height])
    f.Position = [0 0 400, 700];

    imagesc(abs(custom.rss(coil_maps)));
    title('Combined Coil Sensitivity Map', 'FontSize', font_title);
    colorbar();

    % Create the mask
    mask = custom.create_mask(k_space, undersampling_rate);
    mask_image = mask(:, :, 1);

    % Apply the mask to the k-space data
    undersampled_k_space = k_space .* mask;
    undersampled_k_space_image = custom.rss(undersampled_k_space);

    % Convert the undersampled k-space data to images
    undersampled_images = custom.ifft2op(undersampled_k_space);
    undersampled_image = custom.rss(undersampled_images);

    % Use the SENSE algorithm to reconstruct the image
    sense_result = custom.sense(undersampled_images, coil_maps, undersampling_rate);

    % Generate SENSE image using Root Sum of Squares
    sense_image = custom.rss(sense_result);

    % Scalt the SENSE image by the undersampling rate
    sense_image = sense_image .* undersampling_rate;

    original = mat2gray(img_combined);
    undersampled_image = mat2gray(undersampled_image);
    sense_image = mat2gray(sense_image);
    delta = abs(original - sense_image);

    % Create a new figure window
    fig = figure;

    % Set the figure window's size ([left bottom width height])
    fig.Position = [0 0 600 800];

    colormap('gray');
    subplot(2, 2, 1);
    imshow(original);
    title('Original Image', 'FontSize', font_title);
    colorbar();

    axis image;

    subplot(2, 2, 2);
    imshow(undersampled_image);
    title('Undersampled Image', 'FontSize', font_title);
    colorbar();

    axis image;

    subplot(2, 2, 3);
    imshow(sense_image);
    title('SENSE Image', 'FontSize', font_title);
    colorbar();

    axis image;

    subplot(2, 2, 4);
    imshow(delta);
    title('Difference Image', 'FontSize', font_title);
    colorbar();

    axis image;

    text = sprintf('Undersampling Rate: %d, Slice: %d', undersampling_rate, slice);

    sgtitle(text, 'FontSize', 20);

    % Calculate the Mean Squared Error (MSE)
    mse = mean((original(:) - sense_image(:)) .^ 2);

    % Calculate the Peak Signal-to-Noise Ratio (PSNR)
    peakval = max(original(:));
    psnr = 20 * log10(peakval / sqrt(mse));

    % Calculate the Signal-to-Noise Ratio (SNR)
    signal_power = mean(original(:) .^ 2);
    noise_power = mean((original(:) - sense_image(:)) .^ 2);
    snr = 10 * log10(signal_power / noise_power);

    % Display the results
    fprintf('MSE: %.4f\n', mse);
    fprintf('PSNR: %.4f dB\n', psnr);
    fprintf('SNR: %.4f dB\n', snr);

    % Create a string that includes all the results
    result_str = sprintf('MSE: %.4f PSNR: %.4f dB SNR: %.4f dB\n', mse, psnr, snr);

    % Display the results in subtext
    custom.subtext(result_str);

    if is_save
        saveas(fig, fullfile(folder_name, sprintf('fastmri_%s_s_%d_r_%d.png', filename_w, slice, undersampling_rate)));
        close(fig);
    end

    results(result_index).filename = filename;
    results(result_index).slice = slice;
    results(result_index).undersampling_rate = undersampling_rate;
    results(result_index).sense_image = sense_image;
    results(result_index).mse = mse;
    results(result_index).psnr = psnr;
    results(result_index).snr = snr;

    result_index = result_index + 1;

end

% Save the results to a .mat file

if is_save
    save(fullfile(folder_name, sprintf('fastmri_%s_results.mat', filename_w)), 'results');
end
