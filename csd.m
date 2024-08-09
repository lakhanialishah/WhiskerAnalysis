function [iCSD, F_mat_inv, r_LFP] = csd(LFP_mat, param)
% Computes the inverse CSD assuming that the current source density is
% nonzero in a cylindrical region and zero outside, following the set-up in
% Petterson et al 2006: "Current-source density estimation based on
% inversion of electrostatic forward solution: Effects of finite extent of 
% neuronal activity and conductivity discontinuities"
% LFP_mat is a matrix with rows corresponding to recordings on a channel.
% Its dimensions are (#channels) by (#timepoints). 
% param is a struct containing the parameters for iCSD:
%       iCSD_method: delta_source, step_current, offcenter_step_current
%       h: channel separation
%       sigma_c: conductivity (assumed constant, isotropic)
%       R: radius of source cylinder
%       dist_off_center: (offcenter_step_current only) distance of
%       recording electrode to the center of the cylinder. 
% iCSD, F_mat_inv, r_LFP : current sources (iCSD), iCSD =
% -F_mat_inv*LFP_mat; and for the kernel method r_LFP is the reconstructed
% LFP. 
% 
% param = {};
% param.iCSD_method = 'delta_source';
% param.h = 24.2; %um
% param.sigma_c = .3; %what the paper has 
% param.R = 21.65; %um
% param.dist_off_center = 0;
% 
iCSD_method = param.iCSD_method;
% 
r_LFP = [];

if strcmp(iCSD_method, 'kernel1D')
    % Kernel method follows Potworowski (Wójcik) et al 2012
    good_channel_list = param.good_channel_list;
    good_channel_list = reshape(good_channel_list, 1, numel(good_channel_list));
    source_zi = -param.channel_sep*good_channel_list';
    sourceR = param.sourceR;
    xy_r = param.R;
    sigma_c = param.sigma_c;
    lambda = param.lambda;

    %%
    z_mesh = linspace(min(good_channel_list), max(good_channel_list),round(10*length(good_channel_list)))*-param.channel_sep;
    dz = z_mesh(2)-z_mesh(1);


    btilde_i = bsxfun(@(x,y) exp(-(x-y).^2/(2*sourceR^2)), z_mesh, source_zi);
    b_i = 0*btilde_i;
    for ii = 1:size(btilde_i,1)
        for jj = 1:length(z_mesh)
            b_i(ii, jj) = (1/(2*sigma_c))*...
                sum(dz*(sqrt((z_mesh(jj) - z_mesh).^2+xy_r^2) - abs(z_mesh - z_mesh(jj))).*btilde_i(ii,:));
        end
    end

    K_zzprime = b_i*b_i';
    K_ztildez = btilde_i*b_i';


    F_mat_inv = (K_ztildez')*inv(K_zzprime + lambda*eye(size(K_zzprime)));

    iCSD = -F_mat_inv*LFP_mat(good_channel_list, :);

    %% reconstruct LFP using basis functions
    
    r_LFP = -b_i'*iCSD;
    
    inds_good_ch = round(interp1(z_mesh, 1:length(z_mesh), good_channel_list*-param.channel_sep));
    r_LFP = r_LFP(inds_good_ch, :);
    %%
    return;


elseif strcmp(iCSD_method, 'delta_source')
    %%
    h = param.h;
    sigma_c = param.sigma_c;
    R = param.R;
    num_channels = size(LFP_mat, 1);
    
    i_ch = 1:num_channels;
    j_ch = i_ch';
    
    j_minus_i_mat = bsxfun(@minus, -i_ch, -j_ch);
    
    F_mat = (h^2/(2*sigma_c))*(sqrt(j_minus_i_mat.^2 + (R/h)^2) - abs(j_minus_i_mat));
    

    
    
    
elseif strcmp(iCSD_method, 'step_current')

    h = param.h;
    sigma_c = param.sigma_c;
    R = param.R;
    num_channels = size(LFP_mat, 1);
    
    F_mat = zeros(num_channels);
    % because of the symmetry of the integrand, it is not necessary to
    % integrate this for all combinations of i,j. Only unique values
    % of|j-i| require integration
    i_ch = 1:num_channels;
    j_ch = i_ch';
    j_minus_i_mat = bsxfun(@minus, -i_ch, -j_ch);
    j_minus_i_values = unique(abs(j_minus_i_mat));
    for ii_dij = 1:length(j_minus_i_values)
        d_ij = j_minus_i_values(ii_dij);
        F_int_jminusi = integral(@(x) stepFunction_iCSDintegrand(x, R, sigma_c), ...
                h*d_ij-h/2, h*d_ij+h/2);
        F_mat(abs(j_minus_i_mat) == d_ij) = F_int_jminusi;
    end
    
    
elseif strcmp(iCSD_method, 'offcenter_step_current')

    h = param.h;
    sigma_c = param.sigma_c;
    R = param.R;
    num_channels = size(LFP_mat, 1);
    dR = param.dist_off_center;
    
    F_mat = zeros(num_channels);
    % because of the symmetry of the integrand, it is not necessary to
    % integrate this for all combinations of i,j. Only unique values
    % of |j-i| require integration in the z-direction. 
    
    i_ch = 1:num_channels;
    j_ch = i_ch';
    j_minus_i_mat = bsxfun(@minus, -i_ch, -j_ch);
    j_minus_i_values = unique(abs(j_minus_i_mat));
    for ii_dij = 1:length(j_minus_i_values)
        d_ij = j_minus_i_values(ii_dij);
        
        % unfortunately, this integral is 3D with funky limits. Assume
        % (x_e,y_e) of the electrode is at x_e = dR and y_e = 0. 
        ylim_lo = 0;    % symmetry! multiply total integral by 2
        ylim_fun_hi = @(x) sqrt(R^2 - x.^2);
        F_int_jminusi = 2*integral3(@(x, y, z) stepFunctionOffCenter_iCSDintegrand(x-dR, y, z, sigma_c), ...
                -R, R, ylim_lo, ylim_fun_hi, h*d_ij-h/2, h*d_ij+h/2);
        F_mat(abs(j_minus_i_mat) == d_ij) = F_int_jminusi;
    end
    F_mat = F_mat/(2*pi);

end


F_mat_inv = inv(F_mat);
iCSD = F_mat\LFP_mat; %I put in the .