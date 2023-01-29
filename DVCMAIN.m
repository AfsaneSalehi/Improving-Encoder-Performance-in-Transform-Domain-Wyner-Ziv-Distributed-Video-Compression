%--------------------------Codec for Wyner-ziv ---------------------------%
clear;
clc;
addpath(genpath ('.\Videos'));
addpath(genpath ('.\Utilities'));
addpath(genpath ('.\Header'));
addpath(genpath ('.\Residual'));
addpath(genpath ('.'));
addpath('.\Motion Estimation');
addpath('.\Motion Vector Postprocessing');
addpath('.\0 Traditional SI Extrapolation');
addpath('.\1 AR Model');
addpath('.\2 MAAR Model');
addpath('.\Motion Analysis');
addpath('.\Motion Compensation');

% ---------------------Define input video sequence WZ----------------------%
videoName                       = 'foreman_qcif.yuv';
format                          = 'qcif';
width                           = 176;
height                          = 144;
frame_rate                      = 30;
GOP_len                         = 2;
GOP_num                         = 1;
frame_num                       = GOP_len*GOP_num + 1;
init2last                       = [0,frame_num-1];
block_size_WZ                   = 8; %block size
[Y,U,V]                         = yuvRead(videoName, width, height ,frame_num);
level_QP                        = 8;
num_of_bitplanes                = ceil(log2(level_QP));
L                               = 1; %The radius of input feaure 
overlap_size                    = block_size_WZ/2; % The overlapped size of each block

%---------------------------Side Info input-------------------------------%
params.block_size               = 8;
params.search_range             = 16;
params.step_size                = 1;
params.im_rows                  = height;
params.im_cols                  = width;

%---------------------------H.264 input-----------------------------------%
global h w h2 w2  block_size
idx                             = 1;
idx2                            = 1;
block_size                      = 8;

%--------------------Initializes Indication matrix.-----------------------%
indicationMatrix_Y = zeros(height/block_size_WZ,width/block_size_WZ,'double');
indicationMatrix_SI = zeros(height/block_size,width/block_size,'double');

% ---------------------Define input for turbo-----------------------------%
Differntial_DC    = 396;

% Systematic bit arrays
% ---------------------
sys1_y = zeros(num_of_bitplanes,Differntial_DC);
sys2_y = zeros(num_of_bitplanes,Differntial_DC);

% Parity bit arrays
% -----------------
parity1_y = zeros(num_of_bitplanes,Differntial_DC); % use of 1/2 rate constituent encoders assumed
parity2_y = zeros(num_of_bitplanes,Differntial_DC); % use of 1/2 rate constituent encoders assumed

% Interleaver
% -----------
interleaver_y   = zeros(num_of_bitplanes,Differntial_DC);

% De-interleaver 
% --------------
deinterleaver_y = zeros(num_of_bitplanes,Differntial_DC);

% Puncturing (Random puncturing mask)
% -----------------------------------   
puncturing_period=32;
iterations=10;
temp_decoded_bit_plane_y = zeros(Differntial_DC,num_of_bitplanes);
bits_spent_WZ_y_bitplane = zeros(1,num_of_bitplanes);

%-------------------------- H.264 Output----------------------------------%
bitstream_Seq_1                       = '';         % Output bitstream
bitstream_Seq_2                       = '';         % Output bitstream
%----------------------------Reporting------------------------------------%
dec_psnr  = [];
SI_psnr   = [];
bits_WZ   = [];
dec_PreI  = [];
dec_NextI = [];
ssim      = [];
time      = [];
total_bits= 0 ;
total_psnr= 0;
RD_WZ_for_even_frames = [];
total_psnr_even_frames= [];
%---------------------Main Programm---------------------------------------%
% Show sample frames
disp([videoName,': ',num2str(GOP_num),' GOPs with length of ',num2str(GOP_len),' frames'])
figure;
c = 0;  % counter
% Read the video sequence
for iFrame = 1:frame_num
    c = c + 1;
    subplot(6,6,c);
    imshow(Y(:,:,iFrame));
    title(['frame #', num2str(iFrame)]);
end
%% Main Loop.

for i = 2:GOP_len:(frame_num-1)
    
    tic
    %% Encoding H.264 Previouse I Frame
    %----------------------------------------------------------------------
    Frame_ODD_Pre_unit8           = Y(:,:,i-1);
    Frame_ODD_Pre_double          = double(Y(:,:,i-1));
    [h,w,u]                       = size(Frame_ODD_Pre_unit8);
    [bits_h264]                   = header(h,w,level_QP,i-1,i+1 );
    bitstream                     = bits_h264;
    bitstream                     = [bitstream '1111'];
    figure()
    image(Frame_ODD_Pre_unit8)
    colormap(gray(256))
    truesize([2*h 2*w])
    title(['Encoding I Frame',num2str(i-1)])
    pause(0.1)
    
    disp('Encoding I Frame')
    disp('QP=')
    disp(level_QP)
    [Frame_ODD_Pre_double_enc,bits1]        = encode_i_frame(Frame_ODD_Pre_double,level_QP);
    bitstream_Seq_1                         = [bitstream bits1];
    
    %% Encoding H.264 Next I Frame.
    %----------------------------------------------------------------------
	Frame_ODD_Next_unit8     = Y(:,:,i+1);
    Frame_ODD_Next_double    = double(Y(:,:,i+1));
    [h2,w2,u2]               = size(Frame_ODD_Next_unit8);
    figure()
    image(Frame_ODD_Next_unit8)
    colormap(gray(256))
    truesize([2*h2 2*w2])
    title(['Encoding I Frame',num2str(i+1)])
    pause(0.1)
    
    disp('Encoding I Frame')
    disp('QP=')
    disp(level_QP)
    [Frame_ODD_Next_double_enc,bits2]        = encode_i_frame(Frame_ODD_Next_double,level_QP);
    bitstream_Seq_2                         = [bitstream bits2];
    
    save(['bitstream','_enc'],'bitstream')  
	load ./bitstream1_enc.mat; % correct packets
    load ./bitstream2_enc.mat; % correct packets

    disp('Decoding I Frame')

    %% Decoding H.264 Previouse I Frame.
    %----------------------------------------------------------------------
    [h,w,~,~,~,m] = dec_header(bitstream_Seq_1);
    idx = idx + m - 1;

    disp('Decoding I Frame')
    idx = idx + 4;
    [Ceq_r1]=decode_i_frame(idx,bitstream_Seq_1);
    
    Ceq_r_ec1_unit8 = repmat(Ceq_r1, 1, 1, 3);
    Ceq_r1_u8=repmat(Ceq_r1, 1, 1, 3);
    image8Bit_PRE_ODD = uint8(255 * mat2gray(Ceq_r_ec1_unit8));
    Ceq_r1_u8= uint8(255 * mat2gray(Ceq_r1_u8));
    figure()
    pause(1)
    image(image8Bit_PRE_ODD)
    title(['Decoding Frame No. ' num2str(i-1)]);
    colormap(gray(256))
    truesize([2*h 2*w])
    drawnow
    idx=1;
    
    %% Decodin H.264 Next I Frame.
    %----------------------------------------------------------------------        
    [h2,w2,QP,Frame_start,Frame_end,m] = dec_header(bitstream_Seq_2);
    idx2 = idx2 + m - 1;

    disp('Decoding I Frame')
    idx2 = idx2 + 4;
    [Ceq_r2]=decode_i_frame(idx2,bitstream_Seq_2);
    Ceq_r_ec2_unit8 = repmat(Ceq_r2, 1, 1, 3);
    image8Bit_Next_ODD = uint8(255 * mat2gray(Ceq_r_ec2_unit8));
    figure()
    pause(1)
    image(image8Bit_Next_ODD)
    title(['Decoding Frame No. ' num2str(i+1)]);
    colormap(gray(256))
    truesize([2*h2 2*w2])
    drawnow
    idx2=1;
    
    %% SI Interpolation.
    %----------------------------------------------------------------------
    Ref1=Ceq_r1;
    Ref2=Ceq_r2;
    Side_Information_Double=(Ceq_r1+Ceq_r2)/2;
%   Side_Information_Double=double(Side_Information(:,:,i));
    Image_Side_Information = repmat(Side_Information_Double, 1, 1, 1);
    image8Bit_Side_Information = uint8(255 * mat2gray(Image_Side_Information));
    figure();
    imshow(image8Bit_Side_Information);
    
    %% Discrete Cosine Transform (DCT) for SI Frames.
    %----------------------------------------------------------------------
    [rowSI, colSI] = size(Side_Information_Double);
    % DCT works on 8x8 matrix. Thus padding is required if size is not a multiple of 8x8.
    padderSI = zeros(mod(rowSI,8),mod(colSI,8));
    % Padding in case size is not a multiple of 8x8 or 16x16.
    SI_DCTMatrix = [Side_Information_Double; double(padderSI)];
    % Size of padded matrix.
    [rowPaddedSI, colPaddedSI] = size(SI_DCTMatrix);
    
    for m = 0:((rowPaddedSI/8)-1)
       for n = 0:((colPaddedSI/8)-1)
          SI_DCTCurrentFrame((m*8)+1:(m+1)*8, (n*8)+1:(n+1)*8) = dct2(SI_DCTMatrix((m*8)+1:(m+1)*8, (n*8)+1:(n+1)*8));
       end
    end 
    
    %% Uniform Quantization of SI.
    %----------------------------------------------------------------------
    SI_Quantized = Uniform_Q(SI_DCTCurrentFrame,level_QP);
	
    %% Zig-Zag Scan_SI.
    %----------------------------------------------------------------------    
    % Performs Zig-Zag for difference Luma and Chroma components.
    % Passing 2nd argument in zigZag function helps in prevention of 
    % sending those MBs which have MV (0,0) and SAD less than 128.
    
    [currentSI_DC, zigZagCurrentLumaAC_SI] = zigZag(SI_Quantized, indicationMatrix_SI, 1);            
    
    %% Calculates Differential DC coefficients.
    %----------------------------------------------------------------------
    differentialCurrent_SI_DC = differentialCoding(currentSI_DC);
    
    %% coding WZ
    %----------------------------------------------------------------------
    figure();
    imshow(Y(:,:,i));

    
    currentFrameY    = Y(:,:,i); % Extracts Luminance Component. 
    referenceFrameY  = Y(:,:,i-1); % Extracts Luminance Component.
    currentFrameY_Double=double(currentFrameY);   
    title(['Original Current Frame ' int2str(i) ' at ENCODER']);

    %%Intializes Motion Vector Matrix to store the values.
    %----------------------------------------------------------------------     
    motionVectorYValue = zeros(height/block_size_WZ,width/block_size_WZ,'double');
   
    % Initializes Difference Frames for Luma and Chroma components.
    %----------------------------------------------------------------------
    differenceFrameY = zeros(height, width, 'double');
  
    %% Motion Estimation.
    %----------------------------------------------------------------------
    % Extracts Reference and Current Frame Blocks for Motion Estimation.
    %----------------------------------------------------------------------
    for p = 0:((height/8)-1)
        for q = 0:((width/8)-1)

            currentMacroBlock = double(currentFrameY((p*8)+1:(p+1)*8, (q*8)+1:(q+1)*8));

            % Creates Search Window of suitable size.

            if((((q*8)+1)-8)<=0)
                   startColumnSearchWindow = 1;
            else
                   startColumnSearchWindow = (((q*8)+1)-8);
            end

            if((((p*8)+1)-8)<=0)
                   startRowSearchWindow = 1;
            else
                   startRowSearchWindow = (((p*8)+1)-8);
            end

            if(((q+1)*16)+8 >= width)
                   endColumnSearchWindow = width;
            else
                   endColumnSearchWindow = ((q+1)*8)+8;
            end

            if(((p+1)*8)+8 >= height)    
                   endRowSearchWindow = height;
            else
                   endRowSearchWindow = ((p+1)*8)+8;
            end

            referenceSW = double(referenceFrameY(startRowSearchWindow : endRowSearchWindow, startColumnSearchWindow : endColumnSearchWindow));

            % Stores relative X, Y values and difference MB.
            [ xShift, yShift, differenceFrameY((p*8)+1:(p+1)*8, (q*8)+1:(q+1)*8)] = exhaustiveSearchAlgorithm( currentMacroBlock, referenceSW);

            % Calculates Motion Vector.
            motionVectorXValue(p+1,q+1) = (p*8 + 1) - (xShift - 1 + startRowSearchWindow);
            motionVectorYValue(p+1,q+1) = (q*8 + 1) - (yShift - 1 + startColumnSearchWindow);

        end
    end
%     figure;
%     imshow(uint8(differenceFrameY));
%     title(['Y Component of Difference frame for Current Frame No. ' int2str(i) ' at ENCODER']);

    % Plots Motion Vector.
    [x,y] = meshgrid(1:((width/8)),1:((height/8)));
    figure;
    quiver(x,y, motionVectorXValue, motionVectorYValue);
    title(['Motion Vector for Frame No. ' int2str(i)]);

    %% Discrete Cosine Transform (DCT) for WZ Frames.
    %----------------------------------------------------------------------
    [rowY, colY] = size(currentFrameY_Double);
    % DCT works on 8x8 matrix. Thus padding is required if size is not a multiple of 8x8.
    padderY = zeros(mod(rowY,8),mod(colY,8));
    % Padding in case size is not a multiple of 8x8 or 16x16.
    YDCTMatrix = [currentFrameY; double(padderY)];
    % Size of padded matrix.
    [rowPaddedY, colPaddedY] = size(YDCTMatrix);
    % Performs DCT for Luma component, processing 8x8 blocks at a time.
    for m = 0:((rowPaddedY/8)-1)
       for n = 0:((colPaddedY/8)-1)
          YDCTCurrentFrame((m*8)+1:(m+1)*8, (n*8)+1:(n+1)*8) = dct2(YDCTMatrix((m*8)+1:(m+1)*8, (n*8)+1:(n+1)*8));
       end
    end
    
    %% Uniform Quantization.
    %----------------------------------------------------------------------   
    WZQuantized = Uniform_Q(YDCTCurrentFrame,level_QP);
    
    %% Zig-Zag Scan.
    %----------------------------------------------------------------------    
    % Performs Zig-Zag for difference Luma and Chroma components.
    % Passing 2nd argument in zigZag function helps in prevention of 
    % sending those MBs which have MV (0,0) and SAD less than 128.
    [currentY_DC, zigZagCurrentLumaAC] = zigZag(WZQuantized, indicationMatrix_Y, 1);            
    
    %% Calculates Differential DC coefficients.
    %----------------------------------------------------------------------
    differentialCurrent_Y_DC = differentialCoding(currentY_DC);
    
    %% BITPLANE EXTRACTION.
    %----------------------------------------------------------------------
	bitplane_Y = dec2bin(differentialCurrent_Y_DC,num_of_bitplanes);
    bitplane_SI = dec2bin(differentialCurrent_SI_DC,num_of_bitplanes); %FOR SYSTEMATIC BITPLANE
    
	%% turbo encoder.
    %----------------------------------------------------------------------
    for plane = 1:num_of_bitplanes
        data_group_SI = bitplane_SI(:,plane)';
        data_group_SI = str2num(data_group_SI(:))';
        data_group_y = bitplane_Y(:,plane)';
        data_group_y = str2num(data_group_y(:))';
        [parity1_y(plane,:), sys1_y(plane,:), parity2_y(plane,:), sys2_y(plane,:), interleaver_y(plane,:), deinterleaver_y(plane,:)]=...
        turbo_encode(data_group_y);
    
    end
    	
    
    %% Estimation of alpha
    % -------------------
	[alpha, alpha_scaled ]= estimate_alpha_fmo(data_group_y,data_group_SI);
	
    % Find the conditional probabibilties (channel measurement)of 0 and 1
	% for each bitplane of each WZ Frame
	[bitplane_conditional_probability_0_y, bitplane_conditional_probability_1_y] = qpel_bit_probabilities(data_group_SI,alpha,alpha_scaled,level_QP);
	
    %% Decode (turbo) all bitplanes 
    % -----------------------------
    for plane = 1:num_of_bitplanes
            
        current_bitplane_conditional_prob_0_y = bitplane_conditional_probability_0_y(:,:,plane)';
        current_bitplane_conditional_prob_0_y = current_bitplane_conditional_prob_0_y(:)';
        current_bitplane_conditional_prob_1_y = bitplane_conditional_probability_1_y(:,:,plane)';
        current_bitplane_conditional_prob_1_y = current_bitplane_conditional_prob_1_y(:)';

   
        % input to decoder
        par1_y      = parity1_y(plane,:);
        par2_y      = parity2_y(plane,:);
        interl_y    = interleaver_y(plane,:);
        deinterl_y  = deinterleaver_y(plane,:);
        original    = sys1_y(plane,:);
            
        [decoded_current_bit_plane_y,bits_spent_WZ_y_bitplane, BER]=...
                turbo_decode(current_bitplane_conditional_prob_0_y,current_bitplane_conditional_prob_1_y,par1_y, par2_y,... 
                interl_y,deinterl_y,original,iterations, puncturing_period);

    end
    
    %% Inverse Differential Coding.
    %----------------------------------------------------------------------
    decoderCurrentYDC = inverseDifferentialCoding(decoded_current_bit_plane_y);

    %% Inverse Zig-Zag.
    %----------------------------------------------------------------------
    inverseZigZagCurrentY = inverseZigZag(height, width, decoderCurrentYDC, zigZagCurrentLumaAC, indicationMatrix_Y, 1);

    %% Inverse Quantization.
	%----------------------------------------------------------------------    
    y_hatWZ = Uniform_Qinv(inverseZigZagCurrentY,8);
    
    %% Inverse DCT.
	%---------------------------------------------------------------------- 
    [originalRow, originalColumn] = size(y_hatWZ);    
    invDCTLuma = inverseDCT(originalRow, originalColumn, y_hatWZ);  
	
    %% Reconstruction.
	%---------------------------------------------------------------------- 
    [reconstructedMatrix] = reconstruct(invDCTLuma,Side_Information_Double);
    figure; (imshow(uint8(  reconstructedMatrix  )));
    title(['Reconstruct ' int2str(i) ' at ENCODER']); 
    temp1 = toc;
    %% Info for Reporting.
	%---------------------------------------------------------------------- 
    reconstructedMatrix_Double = double(reconstructedMatrix);
    time  = [time,temp1];
    
    temp2 = calculate_psnr(currentFrameY_Double,reconstructedMatrix_Double);
    temp4 = calculate_psnr(currentFrameY_Double,Side_Information_Double);
    temp5 = calculate_psnr (Frame_ODD_Pre_double,Ceq_r1);
    temp6 = calculate_psnr (Frame_ODD_Next_double,Ceq_r2);
    
    
    dec_psnr = [dec_psnr,temp2];
    SI_psnr  = [SI_psnr,temp4];
    dec_PreI = [dec_PreI,temp5];
    dec_NextI = [dec_NextI,temp6];
    
    imwrite(reconstructedMatrix,'.\Transfered\WZFrame.jpg')
    fileInfo                  = dir('.\Transfered\WZFrame.jpg');
    bytesTransfered           = [fileInfo.bytes];
    bitsTransfered            = bytesTransfered*8;
    total_bits                = total_bits+bitsTransfered;
    total_psnr                = total_psnr+temp2;
    RD_WZ_total               = (frame_rate/2)*(total_bits)/1000;
    RD_WZ                     = (frame_rate/2)*(bitsTransfered)/1000;
    bits_WZ                   = [bits_WZ,RD_WZ];    
    RD_WZ_for_even_frames     = [RD_WZ_for_even_frames,RD_WZ_total];
    total_psnr_even_frames    = [total_psnr_even_frames,total_psnr];
    
    
end
%% Evaluation.
%---------------------------------------------------------------------- 
dec_psnr_mean  = mean(dec_psnr);
SI_psnr_mean   = mean(SI_psnr);
RD_WZ_mean     = (frame_rate/2)*mean(bits_WZ)/1000;
time_mean      = mean(time);
disp(['The mean dec Psnr = ',num2str(dec_psnr_mean),' dB']);
disp(['The mean SI Psnr = ',num2str(SI_psnr_mean),' dB']);
disp(['The RD_WZ = ',num2str(RD_WZ_mean),' kbps']);
disp(['The mean Time = ',num2str(time_mean),' s']);
figure;
semilogy(2:2:frame_num,dec_psnr,'k-*','linewidth',2.0);
hold on
semilogy(2:2:frame_num,SI_psnr,'m-o','linewidth',2.0);

title('SI and WZ PSnr');
xlabel('Even Frame Number');ylabel('Psnr');
legend('WZ_psnr','SI_psnr');
grid on
axis tight
hold off

figure;
plot(RD_WZ_for_even_frames,total_psnr_even_frames,'g-o','linewidth',2.0);
title('RD');
xlabel('Rate');ylabel('Psnr');
legend('RD','psnr');
grid on
axis tight
