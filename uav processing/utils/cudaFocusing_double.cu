#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <iostream>

#define checkCudaErrors(val) check_cuda((val), #val, __FILE__, __LINE__)
void check_cuda(cudaError_t result, char const *const func, const char *const file, int const line)
{
    if (result)
    {
        std::cerr << "CUDA error = " << static_cast<unsigned int>(result) << " at " << file << ":" << line << " '" << func << "' \n";
        // Make sure we call CUDA Device Reset before exiting
        cudaDeviceReset();
        mexErrMsgTxt("CUDA Error,exit.");
    }
}

__device__ double gaussActivFunc(double x, double sigma)
{
    double out = exp(-0.5 * (x / sigma) * (x / sigma));
    return out;
}

__device__ double2 linear_interp_comp(double const *x, double2 const *y, double const xq, int const y_zero_idx, int const N)
{
    // RC is the whole matrix, has to take only 1 column
    int const y_end_idx = y_zero_idx + N;

    // Manage extrapolation
    if (xq <= x[0])
        return y[0];
    if (xq >= x[N - 1])
        return y[y_end_idx - 1];
    double2 yout;
    // search the left value
    int j = 0;
    while (x[j + 1] < xq)
        j++;
    double xL = x[j], xR = x[j + 1];
    double2 yL = y[y_zero_idx + j], yR = y[y_zero_idx + j + 1];
    double grad_real = (yR.x - yL.x) / (xR - xL);
    double grad_imag = (yR.y - yL.y) / (xR - xL);
    yout.x = yL.x + grad_real * (xq - xL);
    yout.y = yL.y + grad_imag * (xq - xL);
    return yout;
}

__device__ double2 exp_comp(double const ampl, double const phase)
{
    double2 out;
    out.x = ampl * cos(phase);
    out.y = ampl * sin(phase);
    return out;
}

__device__ double2 mult_comp(double2 const a, double2 const b)
{
    double2 out;
    out.x = a.x * b.x - a.y * b.y;
    out.y = a.x * b.y + a.y * b.x;
    return out;
}
__global__ void focusTDBPKernel(double const *X, double const *Y, double const z0, double const *TX_pos_x,
                                double const *TX_pos_y, double const *TX_pos_z, double const *RX_pos_x, double const *RX_pos_y, double const *RX_pos_z,
                                double const lambda, double const Dk, double2 const *RC, double const *t, double const f0, double const k_rx_0,
                                double2 *Sn, double *Wn, int const N_pixel, int const N_RC, int const tau, int const squint, double const speed_weight)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= N_pixel)
        return;
    int out_i = i + N_pixel*squint;
    
    double X_i = X[i];
    double RX_pos_x_tau = RX_pos_x[tau];
    if (X_i < RX_pos_x_tau)
    {
        // Backlobe of antenna, pixel is 0
        Wn[out_i] += 0.0;
        Sn[out_i].x += 0.0;
        Sn[out_i].y += 0.0;
        return;
    };
    double const C = 2.99792458e8;
    double const pi = 3.1415926535897932385;
    double Y_i = Y[i];
    double TX_pos_x_tau = TX_pos_x[tau];
    double TX_pos_y_tau = TX_pos_y[tau];
    double TX_pos_z_tau = TX_pos_z[tau];
    double RX_pos_y_tau = RX_pos_y[tau];
    double RX_pos_z_tau = RX_pos_z[tau];

    // Range distances from the tx antenna [m]
    double R_tx = sqrt((TX_pos_x_tau - X_i) * (TX_pos_x_tau - X_i) + (TX_pos_y_tau - Y_i) * (TX_pos_y_tau - Y_i) + (TX_pos_z_tau - z0) * (TX_pos_z_tau - z0));
    // Range distances from the rx antenna [m]
    double R_rx = sqrt((RX_pos_x_tau - X_i) * (RX_pos_x_tau - X_i) + (RX_pos_y_tau - Y_i) * (RX_pos_y_tau - Y_i) + (RX_pos_z_tau - z0) * (RX_pos_z_tau - z0));
    // Total Tx-target-Rx distance [m]
    double distance = R_tx + R_rx;
    double delay = distance / C;

    // Compute target wave number
    double R = sqrt((RX_pos_x_tau - X_i) * (RX_pos_x_tau - X_i) + (RX_pos_y_tau - Y_i) * (RX_pos_y_tau - Y_i));
    double psi = asin((Y_i - RX_pos_y_tau) / R);
    double k_rx = sin(psi) * 2 * pi / lambda;

    // Weight function
    double sigma = Dk / 2;
    
    double Wn_i = speed_weight * gaussActivFunc(k_rx - k_rx_0, sigma);

    // Backprojection of data from a single Radar position
    int const RC_zero_idx = tau * N_RC;
    double2 RC_1 = linear_interp_comp(t, RC, delay, RC_zero_idx, N_RC);

    double2 RC_2 = mult_comp(RC_1, exp_comp(1, 2 * pi * f0 * delay));
    Sn[out_i].x += Wn_i * RC_2.x;
    Sn[out_i].y += Wn_i * RC_2.y;
    Wn[out_i] += Wn_i;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //========================================= Input variables
    mxGPUArray const *X, *Y, *RC, *t, *TX_pos_x, *TX_pos_y, *TX_pos_z, *RX_pos_x, *RX_pos_y, *RX_pos_z;

    double const *d_X, *d_Y, *d_t, *d_TX_pos_x, *d_TX_pos_y, *d_TX_pos_z, *d_RX_pos_x, *d_RX_pos_y, *d_RX_pos_z,*k_rx_0_vec,*RX_speed;
    double2 const *d_RC;
    //========================================= Constants
    double const *_z0, *_lambda, *_Dk, *_f0,*_median_speed;

    //========================================= Output variables
    mxGPUArray *Sn, *Wn;
    double2 *d_Sn;
    double *d_Wn;
    
    int N_pixel, N_RC, N_tau, N_squints;
    char const *const errId = "parallel:gpu:mexGPUExample:InvalidInput";

    //========================================= Input validation
    if (nrhs != 17 || nlhs != 2)
        mexErrMsgIdAndTxt(errId, "Wrong number of input/output arguments.");
    if (!(mxIsGPUArray(prhs[0]) && mxIsGPUArray(prhs[1])))
        mexErrMsgIdAndTxt(errId, "Input must be GPUArray");

    int const threadsPerBlock = 256;
    int blocksPerGrid;

    mxInitGPU();
    //========================================= Initialize variables
    X = mxGPUCreateFromMxArray(prhs[0]);
    Y = mxGPUCreateFromMxArray(prhs[1]);
    _z0 = (double const *)mxGetData(prhs[2]);
    TX_pos_x = mxGPUCreateFromMxArray(prhs[3]);
    TX_pos_y = mxGPUCreateFromMxArray(prhs[4]);
    TX_pos_z = mxGPUCreateFromMxArray(prhs[5]);
    RX_pos_x = mxGPUCreateFromMxArray(prhs[6]);
    RX_pos_y = mxGPUCreateFromMxArray(prhs[7]);
    RX_pos_z = mxGPUCreateFromMxArray(prhs[8]);
    _lambda = (double const *)mxGetData(prhs[9]);
    _Dk = (double const *)mxGetData(prhs[10]);
    RC = mxGPUCreateFromMxArray(prhs[11]);
    t = mxGPUCreateFromMxArray(prhs[12]);
    _f0 = (double const *)mxGetData(prhs[13]);
    k_rx_0_vec = (double const *) mxGetData(prhs[14]);
    RX_speed = (double const *) mxGetData(prhs[15]);
    _median_speed = (double const *) mxGetData(prhs[16]);

    double const z0 = _z0[0];
    double const lambda = _lambda[0];
    double const Dk = _Dk[0];
    double const f0 = _f0[0];
    double const median_speed = _median_speed[0];

    if (mxGPUGetClassID(X) != mxDOUBLE_CLASS)
    {
        mexErrMsgIdAndTxt(errId, "Input must be double");
    }

    //========================================= Initialize pointers
    d_X = (double const *)mxGPUGetDataReadOnly(X);
    d_Y = (double const *)mxGPUGetDataReadOnly(Y);
    d_TX_pos_x = (double const *)mxGPUGetDataReadOnly(TX_pos_x);
    d_TX_pos_y = (double const *)mxGPUGetDataReadOnly(TX_pos_y);
    d_TX_pos_z = (double const *)mxGPUGetDataReadOnly(TX_pos_z);
    d_RX_pos_x = (double const *)mxGPUGetDataReadOnly(RX_pos_x);
    d_RX_pos_y = (double const *)mxGPUGetDataReadOnly(RX_pos_y);
    d_RX_pos_z = (double const *)mxGPUGetDataReadOnly(RX_pos_z);
    d_RC = (double2 const *)mxGPUGetDataReadOnly(RC);
    d_t = (double const *)mxGPUGetDataReadOnly(t);

    //========================================= Create ouput array
    N_pixel = (int)mxGPUGetNumberOfElements(X);
    N_RC = (int)mxGPUGetDimensions(RC)[0];
    N_tau = (int)mxGPUGetDimensions(RC)[1];
    N_squints = (int)mxGetDimensions(prhs[14])[0];
    
    mwSize out_N_dim = 3, out_dims[3] = {mxGPUGetDimensions(X)[0],mxGPUGetDimensions(X)[1],N_squints}; 
    
    Sn = mxGPUCreateGPUArray(out_N_dim,
                             out_dims,
                             mxGPUGetClassID(RC),
                             mxGPUGetComplexity(RC),
                             MX_GPU_INITIALIZE_VALUES);
    d_Sn = (double2 *)mxGPUGetData(Sn);
    Wn = mxGPUCreateGPUArray(out_N_dim,
                             out_dims,
                             mxGPUGetClassID(X),
                             mxGPUGetComplexity(X),
                             MX_GPU_INITIALIZE_VALUES);
    d_Wn = (double *)mxGPUGetData(Wn);

    //========================================= Elaboration

    
    std::cout << "\nN_pix " << N_pixel << ", N_RC " << N_RC << "\n";
    std::cout << "lambda " << lambda << " freq " << f0 << "\n";
    std::cout << "RC dim 1 " << N_RC << "\n";
    std::cout << "RC dim 2 " << N_tau << "\n";
    std::cout << "N squints " << N_squints << "\n";
    
    blocksPerGrid = (N_pixel + threadsPerBlock - 1) / threadsPerBlock;
    checkCudaErrors(cudaPeekAtLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    double speed_weight;
    for(int squint = 0; squint < N_squints;squint++){
        for (int tau = 0; tau < N_tau; tau++)
        {
            speed_weight = RX_speed[tau] / median_speed;
            focusTDBPKernel<<<blocksPerGrid, threadsPerBlock>>>(
                d_X, d_Y, z0, d_TX_pos_x, d_TX_pos_y, d_TX_pos_z, d_RX_pos_x, d_RX_pos_y, d_RX_pos_z,
                lambda, Dk, d_RC, d_t, f0, k_rx_0_vec[squint], d_Sn, d_Wn, N_pixel, N_RC, tau, squint, speed_weight);
        }
        std::cout << "Squint k_rx " << k_rx_0_vec[squint] << "\n";
        std::cout << "Squint n " << squint +1 << " of " << N_squints << "\n";
    }
    checkCudaErrors(cudaPeekAtLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    // OUTPUT

    plhs[0] = mxGPUCreateMxArrayOnGPU(Sn);
    plhs[1] = mxGPUCreateMxArrayOnGPU(Wn);

    // Destroy
    mxGPUDestroyGPUArray(Sn);
    mxGPUDestroyGPUArray(Wn);
    mxGPUDestroyGPUArray(X);
    mxGPUDestroyGPUArray(Y);
    mxGPUDestroyGPUArray(TX_pos_x);
    mxGPUDestroyGPUArray(TX_pos_y);
    mxGPUDestroyGPUArray(TX_pos_z);
    mxGPUDestroyGPUArray(RX_pos_x);
    mxGPUDestroyGPUArray(RX_pos_y);
    mxGPUDestroyGPUArray(RX_pos_z);
    mxGPUDestroyGPUArray(RC);
    mxGPUDestroyGPUArray(t);
    //mxGPUDestroyGPUArray(k_rx_0_vec);
}