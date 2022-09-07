#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <iostream>
#include <string.h>
#include "uEye.h"
/* Function call syntax
 *  Initialize camera: uEye_camera(0) 
 *  Take a single frame image: I = uEye_camera(1)  
 *  Start video: uEye_camera(2)
 *  Fetch a single frame from stored video images: I = uEye_camera(3) 
 *  Close camera: uEye_camera(4)
 *  Get frame rate: R = uEye_camera(5)
 *  Set frame rate: uEye_camera(6, R)
 *  Get exposure times: uEye_camera(7)
 *  Set exposure times: uEye_camera(8, T)
 */
struct cam_props {
    HIDS hCam;
    INT nSizeX;
    INT nSizeY;
    INT nBitsPerPixel;
    INT		lMemId;			// camera main memory - buffer ID
    char*	pcImgMem;			// camera main memory - pointer to buffer
    INT nSizeXMem;
    INT nSizeYMem;
    INT nBitsPerPixelMem;
    INT nPitchMem;
    INT nBytesPerPixel;
    INT nImageSize;
};

struct cam_props cp;
bool isopen = false;
INT timeout = 50; // 500 ms
void convToMxArray(const char * frame, mxArray *M);
void GetBitsPerPixel(HIDS hCam, INT * pnBitsPerPixel);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int mode = mxGetScalar(prhs[0]);
    int err = 0;
    INT nRet;
    switch (mode){
        case 0: // Initialize camera
        {
            INT nNumCam;
            is_GetNumberOfCameras(&nNumCam);
            if (nNumCam == 0) {
                std::cout << "No camera found.\n";
                err =1;
                plhs[0] = mxCreateDoubleScalar(err);
                return;
            } else {
                std::cout << nNumCam << " camera(s) found.\n";
            }
            HIDS hCam = 0;
            nRet = is_InitCamera (&hCam, NULL);
            if (nRet != IS_SUCCESS){
                std::cout << "Initializing camera NOT successful.\n";
                //Check if GigE uEye SE needs a new starter firmware
                if (nRet == IS_STARTER_FW_UPLOAD_NEEDED){
                    //Calculate time needed for updating the starter firmware
                    INT nTime;
                    is_GetDuration (hCam, IS_SE_STARTER_FW_UPLOAD, &nTime);
                    std::cout << "Firmware being updated.\n";
                    //Upload new starter firmware during initialization
                    hCam =  hCam | IS_ALLOW_STARTER_FW_UPLOAD;
                    nRet = is_InitCamera (&hCam, NULL);
                }
                else {
                    std::cout << "Initializing camera failed!\n" << "Return code = " << nRet << std::endl;
                    err = 1;
                    plhs[0] = mxCreateDoubleScalar(err);
                    return;
                }
            } else {
                std::cout << "Initializing camera successful.\n";
            }
            SENSORINFO SensorInfo;
            INT nSizeX, nSizeY;
            nRet = is_GetSensorInfo(hCam, &SensorInfo);
            if (nRet == IS_SUCCESS){
                std::cout << "Sensor color mode: " << int(SensorInfo.nColorMode) << std::endl;
                nSizeX = SensorInfo.nMaxWidth;
                nSizeY = SensorInfo.nMaxHeight;
                std::cout << "Max width: " << nSizeX << std::endl;
                std::cout << "Max height: " << nSizeY << std::endl;
            } else {
                std::cout << "!! Getting sensor info failed !!\n";
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
                return;
            }
            INT nBitsPerPixel;
            GetBitsPerPixel(hCam, &nBitsPerPixel);
            std::cout << "Number of bits per pixel = " << nBitsPerPixel << std::endl;
            // allocate memory
            INT		lMemId;			// camera main memory - buffer ID
            char*	pcImgMem;			// camera main memory - pointer to buffer
            nRet = is_AllocImageMem(hCam, nSizeX, nSizeY, nBitsPerPixel, &pcImgMem, &lMemId);
            if (nRet == IS_SUCCESS){
                std::cout << "Allocating memory successful.\n";
            } else {
                std::cout << "!! Allocating memory failed !!\n";
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
                return;
            }
            is_SetImageMem(hCam, pcImgMem, lMemId ); // set memory active
            INT nSizeXMem, nSizeYMem, nBitsPerPixelMem, nPitchMem;
            is_InquireImageMem (hCam, pcImgMem, lMemId, &nSizeXMem, &nSizeYMem, &nBitsPerPixelMem, &nPitchMem ); // read back for calculation
            
            INT nBytesPerPixel = (nBitsPerPixelMem+1)/8; // calculate bytes per pixel
            INT nImageSize = nSizeXMem * nSizeYMem * nBytesPerPixel; // image size [bytes]
            
            std::cout << "Reading back memory properties\n" << "Size X: " << nSizeXMem << std::endl
            << "Size Y: " << nSizeYMem << std::endl << "Bits per pixel: " << nBitsPerPixelMem << std::endl
            << "Pitch: " << nPitchMem << std::endl << "Bytes per pixel: " << nBytesPerPixel << std::endl
            << "Image size: " << nImageSize << std::endl;
             is_SetDisplayMode(hCam, IS_SET_DM_DIB);
            // store in cp
            cp.hCam = hCam;
            cp.nSizeX = nSizeX;
            cp.nSizeY = nSizeY;
            cp.nBitsPerPixel = nBitsPerPixel;
            cp.lMemId =  lMemId;
            cp.pcImgMem = pcImgMem;
            cp.nSizeXMem = nSizeXMem;
            cp.nSizeYMem = nSizeYMem;
            cp.nBitsPerPixelMem = nBitsPerPixelMem;
            cp.nPitchMem = nPitchMem;
            cp.nBytesPerPixel = nBytesPerPixel;
            cp.nImageSize = nImageSize;
            
            plhs[0] = mxCreateDoubleScalar(err);
            plhs[1] = mxCreateDoubleScalar(nSizeXMem);
            plhs[2] = mxCreateDoubleScalar(nSizeYMem);
            isopen = true;
            break;
        }
        case 1:
        { // take a single frame
            if (!isopen){
                std::cout << "Camera is not open. \n";
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
                return;
            }
            
            nRet = is_FreezeVideo(cp.hCam, timeout);
            if (nRet != IS_SUCCESS) {
                std::cout << "!! Image capture failed !!\n";
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
                return;
            }
            mxArray * frameData = mxCreateDoubleMatrix(1, cp.nImageSize, mxREAL);
            convToMxArray(cp.pcImgMem, frameData);
            plhs[0] = mxCreateDoubleScalar(err);
            plhs[1] = frameData;
            plhs[2] = mxCreateDoubleScalar(cp.nSizeXMem);
            plhs[3] = mxCreateDoubleScalar(cp.nSizeYMem);
            break;
        }
        case 2:
        { // start video
            if (!isopen){
                std::cout << "Camera is not open.\n";
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
                return;
            }
            // start live video
            nRet = is_CaptureVideo(cp.hCam, IS_WAIT);
            std::cout << "Starting video.\n";
            plhs[0] = mxCreateDoubleScalar(err);
            break;
        }
        case 3:
        { // fetch image data
            if (!isopen){
                std::cout << "!! Camera is not open !!\n";
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
                return;
            }
            mxArray * frameData = mxCreateDoubleMatrix(1, cp.nImageSize, mxREAL);
            convToMxArray(cp.pcImgMem, frameData);
            plhs[0] = mxCreateDoubleScalar(err);
            plhs[1] = frameData;
            plhs[2] = mxCreateDoubleScalar(cp.nSizeXMem);
            plhs[3] = mxCreateDoubleScalar(cp.nSizeYMem);
            break;
        }
        case 4:
        { // close camera
            is_FreeImageMem(cp.hCam, cp.pcImgMem, cp.lMemId);
            nRet = is_ExitCamera(cp.hCam);
            if (nRet != IS_SUCCESS){
                std::cout << "Exiting camera failed!\n" << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                std::cout << "Exiting camera.\n";
                plhs[0] = mxCreateDoubleScalar(err);
            }
            break;
        }
        case 5:
        { // get frame rate
            double FR;
            nRet = is_GetFramesPerSecond (cp.hCam, &FR);
            if (nRet != IS_SUCCESS){
                std::cout << "Getting frame rate failed!\n" << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                std::cout << "Current frame rate = " << FR << std::endl;
                plhs[0] = mxCreateDoubleScalar(FR);
            }
            break;
        }
        case 6:
        { // set frame rate
            double FR = mxGetScalar(prhs[1]);
            double newFR;
            nRet = is_SetFrameRate(cp.hCam, FR, &newFR);
            if (nRet != IS_SUCCESS){
                std::cout << "Setting frame rate failed!\n" << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                std::cout << "Current frame rate = " << newFR << std::endl;
                plhs[0] = mxCreateDoubleScalar(newFR);
            }
            break;
        }
        case 7:
        { // get exposure time
            double exp_time;
            double exp_time_max;
            double exp_time_min;
            double exp_time_inc;
            // current setting
            nRet = is_Exposure(cp.hCam, IS_EXPOSURE_CMD_GET_EXPOSURE, &exp_time, sizeof(double));
            if (nRet != IS_SUCCESS){
                std::cout << "Getting exposure time failed!\n" << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                plhs[0] = mxCreateDoubleScalar(exp_time);
                std::cout << "Current exposure time = " << exp_time <<std::endl;
            }
            // max
            nRet = is_Exposure(cp.hCam, IS_EXPOSURE_CMD_GET_EXPOSURE_RANGE_MAX, &exp_time_max, sizeof(double));
            if (nRet != IS_SUCCESS){
                std::cout << "Getting exposure time failed!\n" << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                plhs[1] = mxCreateDoubleScalar(exp_time_max);
                std::cout << "Max exposure time = " << exp_time_max <<std::endl;
            }
            // min
            nRet = is_Exposure(cp.hCam, IS_EXPOSURE_CMD_GET_EXPOSURE_RANGE_MIN, &exp_time_min, sizeof(double));
            if (nRet != IS_SUCCESS){
                std::cout << "Getting exposure time failed!\n" << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                plhs[2] = mxCreateDoubleScalar(exp_time_min);
                std::cout << "Min exposure time = " << exp_time_min <<std::endl;
            }
            // increment
            nRet = is_Exposure(cp.hCam, IS_EXPOSURE_CMD_GET_EXPOSURE_RANGE_INC, &exp_time_inc, sizeof(double));
            if (nRet != IS_SUCCESS){
                std::cout << "Getting exposure time failed!\n" << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                plhs[3] = mxCreateDoubleScalar(exp_time_inc);
                std::cout << "Exposure time increment = " << exp_time_inc <<std::endl;
            }
            break;
        }
        case 8:
        { // set exposure time
            double exp_time = mxGetScalar(prhs[1]);
            std::cout << "Set value = " << exp_time << std::endl;
            nRet = is_Exposure(cp.hCam, IS_EXPOSURE_CMD_SET_EXPOSURE, &exp_time, sizeof(double));
            if (nRet != IS_SUCCESS){
                std::cout << "Setting exposure time failed!\n" << "Set value = " << exp_time << std::endl << "Return code = " << nRet << std::endl;
                err = 1;
                plhs[0] = mxCreateDoubleScalar(err);
            } else {
                plhs[0] = mxCreateDoubleScalar(exp_time);
                std::cout << "Current exposure time =" << exp_time <<std::endl;
            }
            break;
        }
    }
    return;
}

void GetBitsPerPixel(HIDS hCam, INT * pnBitsPerPixel)
{
    switch(is_SetColorMode(hCam, IS_GET_COLOR_MODE))
        {
        case IS_CM_RGBA12_UNPACKED:
        case IS_CM_BGRA12_UNPACKED:
            *pnBitsPerPixel = 64;
            break;

        case IS_CM_RGB12_UNPACKED:
        case IS_CM_BGR12_UNPACKED:
        case IS_CM_RGB10_UNPACKED:
        case IS_CM_BGR10_UNPACKED:
            *pnBitsPerPixel = 48;
            break;

        case IS_CM_RGBA8_PACKED:
        case IS_CM_BGRA8_PACKED:
        case IS_CM_RGB10_PACKED:
        case IS_CM_BGR10_PACKED:
        case IS_CM_RGBY8_PACKED:
        case IS_CM_BGRY8_PACKED:
            *pnBitsPerPixel = 32;
            break;

        case IS_CM_RGB8_PACKED:
        case IS_CM_BGR8_PACKED:
            *pnBitsPerPixel = 24;
            break;

        case IS_CM_BGR565_PACKED:
        case IS_CM_UYVY_PACKED:
        case IS_CM_CBYCRY_PACKED:
            *pnBitsPerPixel = 16;
            break;

        case IS_CM_BGR5_PACKED:
            *pnBitsPerPixel = 15;
            break;

        case IS_CM_MONO16:
        case IS_CM_SENSOR_RAW16:
        case IS_CM_MONO12:
        case IS_CM_SENSOR_RAW12:
        case IS_CM_MONO10:
        case IS_CM_SENSOR_RAW10:
            *pnBitsPerPixel = 16;
            break;

        case IS_CM_RGB8_PLANAR:
            *pnBitsPerPixel = 24;
            break;
        
        case IS_CM_MONO8:
        case IS_CM_SENSOR_RAW8:
        default:
            *pnBitsPerPixel = 8;
            break;
        }
}

void convToMxArray(const char * frame, mxArray *M)
{
    const mwSize *dim = mxGetDimensions(M);
    int nrows = int(dim[0]), ncols = int(dim[1]);
    int numele = ncols*nrows;
    if (ncols != 1 && nrows !=1) mexErrMsgTxt("Invalid matrix size.");
    double *ptrM = mxGetPr(M); 
    for (int n=0;n < numele; n++ ) ptrM[n] = double(frame[n]);
}