#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <memory>
#include <cstdlib>
#include <utility>
#include <cstring>
#include <complex>
#include "utils_verbose_support.h"
#include "utils_time_measure.h"
#include "utils_msg_logger.h"
#include "vie_operator_factory.h"
#include "vie_types.h"
#include "vie_operator_params.h"
#include "geo_spatial_domain.h"
#include "mlaa_base.h"
#include "mex.h"

using  std::cout;
using  std::endl;
using  std::unique_ptr;
using  std::setw;
using  std::setprecision;
using  std::string;
using  std::real;
using  std::complex;

using  marie_core_0_1::Verbose;
using  marie_core_0_1::OperatorFactory;
using  marie_core_0_1::GreenFunction;
using  marie_core_0_1::OperatorStructure;
using  marie_core_0_1::SingularMethod;

using  marie_core_0_1::OperatorParameters;
using  marie_core_0_1::boost_cmplx_4d;
using  marie_core_0_1::axis_size;
using  marie_core_0_1::msg_logger;
using  marie_core_0_1::MessageStreamMode;
using  marie_core_0_1::nothingIdObject;
using  marie_core_0_1::TmFormatType;
using  marie_core_0_1::TimeMeasure;
using  marie_core_0_1::boost_dbl_4d;
using  marie_core_0_1::SpatialDomain;
using  marie_core_0_1::range;
using  marie_core_0_1::const_view_dbl_4to2d;
using  marie_core_0_1::const_view_dbl_4to1d;
using  marie_core_0_1::range;


// functions look for min and max values
int getMinMax_X(const boost_dbl_4d& r, const int* dim, double& min_x, double& max_x);
int getMinMax_Y(const boost_dbl_4d& r, const int* dim, double& min_y, double& max_y);
int getMinMax_Z(const boost_dbl_4d& r, const int* dim, double& min_z, double& max_z);




void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	// Global settings for (debug) output. Must be set before msg_logger call.
    // In fact if it is not called - just default settings are used,
    // See utils_verbose_support.h/cpp  for details
    Verbose::DebugMemory::setShowFlag(false);
    Verbose::UserSupport::setShowFlag(true);
    Verbose::TechnicalDetails::setShowFlag(true);
    Verbose::TimeStamps::setShowFlag(true);
	//    Verbose::setSilentMode();

    // Time format for output
    TimeMeasure::set_default_stamp_format(TmFormatType::yyyymmdd_hhmmss);
    TimeMeasure::set_default_ssdiff_format(TmFormatType::mmssms);
    TimeMeasure::set_time_separator(":");

    
	
    const int number_of_dims = (int)mxGetNumberOfDimensions(prhs[0]);
    const int *dim_array = (int*)mxGetDimensions(prhs[0]);
    const int nelem = (int)mxGetNumberOfElements(prhs[0]);

    double* data = (double*) mxGetPr(prhs[0]);
    double freq = mxGetScalar(prhs[1]);
    int M = mxGetM(prhs[2]);
    int N = mxGetN(prhs[2]);
    char operator_type = * mxArrayToString(prhs[2]);

    char type;
    string singular;

    switch (nrhs)
    {
    case 4:
        type = * mxArrayToString(prhs[3]);
        singular = "DIRECTFN";
        break;
    case 5:
        type = * mxArrayToString(prhs[3]);
        singular = string(mxArrayToString(prhs[4]));
        break;
    default:
        type = 'F';
        singular = "DIRECTFN";
        break;
    }
    cout<<"frequency: "<<freq<<endl;

    if (M * N > 1 )
    {
        const string err_str(" getOPERATORS error: wrong type of operator! ");
        msg_logger().print_message(nothingIdObject(), err_str);
        return;
    }

    double res = 0.0;

    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;

    // Read body
    int count = 0;

    boost_dbl_4d body(boost::extents[dim_array[1]][dim_array[2]][dim_array[3]][dim_array[0]], boost::fortran_storage_order());

    for (int l =0; l<dim_array[3]; ++l)
        for (int k = 0; k<dim_array[2]; ++k)
            for (int j = 0; j < dim_array[1]; ++j)
                for(int i = 0; i< dim_array[0]; ++i)
                {
                    body[i][j][k][l] = data[count];
                    count++;
                }

    res = body[1][0][0][0] - body[0][0][0][0];
    cout<<"Res: "<<res<<endl;

    getMinMax_X(body, dim_array, x_min, x_max);
    getMinMax_Y(body, dim_array, y_min, y_max);
    getMinMax_Z(body, dim_array, z_min, z_max);


	// Message logger setup
    msg_logger().set_output_mode(MessageStreamMode::toConsole);
    msg_logger().print_message(nothingIdObject(), std::string("\n*** VIE Assembly test starts... ***")); 

    unique_ptr<OperatorFactory> upOprFactory(new OperatorFactory("OpFactory", MessageStreamMode::toConsole));

    switch (operator_type)
    {
    case 'K':
        upOprFactory->set_operator_type(GreenFunction::K_type);
        break;

    case 'N':
        upOprFactory->set_operator_type(GreenFunction::N_type);
        break;

    default:
        const string err_str(" getOPERATORS error: wrong type of operator! ");
        msg_logger().print_message(nothingIdObject(), err_str);
        return;
    }
    //    upOprFactory->set_operator_structure(OperatorStructure::Toeplitz);
    //    upOprFactory->set_operator_structure(OperatorStructure::Circulant);

    switch(type)
    {
    case 'T':

        upOprFactory->set_operator_structure(OperatorStructure::Toeplitz);
        cout<<"Toeplitz\n";
        break;

    case 'C':

        upOprFactory->set_operator_structure(OperatorStructure::Circulant);
        cout<<"Circulant\n";
        break;

    default:
        upOprFactory->set_operator_structure(OperatorStructure::FFT_Circulant);
        break;
    }

    switch (strcmp(singular.c_str(),"DEMCEM"))
    {
    	case 0:
    		upOprFactory->set_singular_method(SingularMethod::Demcem);
    		break;
    	default:
    		upOprFactory->set_singular_method(SingularMethod::Directfn);
    		break;
    }
    
    upOprFactory->set_working_frequency(freq);

    unique_ptr<SpatialDomain> upR(new SpatialDomain(string("SpatialDomain"), MessageStreamMode::toConsole));

    upR->set_x_limits(x_min, x_max); //0.5);
    upR->set_y_limits(y_min, y_max); //0.5);
    upR->set_z_limits(z_min, z_max); //0.5);
    upR->set_resolution(res);
    upR->generate();
    upOprFactory->set_domain(upR.get());

    // Calculation Parameters (far, medium, nearby)
    unique_ptr<OperatorParameters> up_op_prm(new OperatorParameters());
    up_op_prm->set_Np_quad_1D_far(4);
    up_op_prm->set_Np_quad_1D_medium(10);
    up_op_prm->set_Np_quad_1D_near(20);
    up_op_prm->set_Np_Directfn(20);
    up_op_prm->set_n_loop_medium(5);
    up_op_prm->set_n_loop_nearby(2);
    upOprFactory->set_op_params(up_op_prm.get());

    // Return a 4d-boost_marray for K or N operator.
    unique_ptr<boost_cmplx_4d> up_bma_c4d(upOprFactory->get_operator(MessageStreamMode::toConsole));
    // Check the proper initialization
    if (!up_bma_c4d) {
        const string err_str(" VIE assembly Error! The Complex-4D operator has not been initialized.");
        msg_logger().print_message(nothingIdObject(), err_str);
        return;
    }
    // Success
    msg_logger().print_message(nothingIdObject(), std::string("\n*** VIE Assembly test successfully passed. ***"));
	
    count = 0;

    // C- storage order()
    boost_cmplx_4d & tmp_c4d_ref = *up_bma_c4d.get();

    int dim_0 = tmp_c4d_ref.shape()[0];
    int dim_1 = tmp_c4d_ref.shape()[1];
    int dim_2 = tmp_c4d_ref.shape()[2];
    int dim_3 = tmp_c4d_ref.shape()[3];

    const int dims[4] = {dim_1, dim_2, dim_3, dim_0};

    plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxCOMPLEX);

    double* op_real = mxGetPr(plhs[0]);
    double* op_imag = mxGetPi(plhs[0]);

    // Read operator from c_storage_order to fortran_storage_order
    // send results to output

    for(int i = 0; i < dim_0; ++i)
    {
        for(int l = 0; l < dim_3; ++l)
        {
            for(int k = 0; k < dim_2; ++k)
            {
                for(int j = 0; j < dim_1; ++j)
                {
                    op_real[count] = real(tmp_c4d_ref[i][j][k][l]);
                    op_imag[count] = imag(tmp_c4d_ref[i][j][k][l]);
                    count++;

                }
            }
        }
    }

}


//------------------------------------------------------------------------------------
// get X-limits of the domain

int getMinMax_X(const boost_dbl_4d& r, const int* dim, double& min_x, double& max_x)
{
	int dim_x = dim[0];
	const_view_dbl_4to1d x_view = r[boost::indices[range(0,dim_x)][0][0][0]];

	max_x = x_view[0];
    min_x = x_view[0];

	for (int i=0; i<dim_x; ++i)
	{	
		if (max_x < x_view[i])
			max_x = x_view[i];

		if(min_x > x_view[i])
			min_x = x_view[i];
	}

	//cout
	
	return 0;
}

//------------------------------------------------------------------------------------
// get Y-limits of the domain

int getMinMax_Y(const boost_dbl_4d& r, const int* dim, double& min_y, double& max_y)
{	
	int dim_y = dim[1];
	const_view_dbl_4to1d y_view = r[boost::indices[0][range(0,dim_y)][0][1]];

	max_y = y_view[0];
    min_y = y_view[0];

	for (int i=0; i<dim_y; ++i)
	{	
		if (max_y < y_view[i])
			max_y = y_view[i];

		if(min_y > y_view[i])
			min_y = y_view[i];
	}
	
	return 0;
}


//------------------------------------------------------------------------------------
// get Z-limits of the domain
int getMinMax_Z(const boost_dbl_4d& r, const int* dim, double& min_z, double& max_z)
{	
	int dim_z = dim[2];

	const_view_dbl_4to1d z_view = r[boost::indices[0][0][range(0,dim_z)][2]];

	max_z = z_view[0];
    min_z = z_view[0];

	for (int i=0; i<dim_z; ++i)
	{	
		if (max_z < z_view[i])
			max_z = z_view[i];

		if(min_z > z_view[i])
			min_z = z_view[i];
	}
	
	return 0;
}
