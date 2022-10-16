#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#include "catch.hpp"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Prob.h"

#include <numeric>
#include <cmath>

// SET THE DEFUALT TOLERANCE FOR TESTS
const float tolerance = 1e-8;

/**
 * This function calculates the l2 norm of the difference between two
 * MultiResol instances.
 * @param[in] MultiResol object 1.
 * @param[in] MultiResol object 2.
 * @return True or false.
 */
bool test_diff(MultiResol &obj1, MultiResol &obj2){

  std::vector<float> coeff_diff;

  for (int i=0; i < obj1.nbr_scale(); i++)
    for (int j=0; j < obj1.size_ima_nl(); j++)
      for (int k=0; k < obj1.size_ima_nc(); k++)
        coeff_diff.push_back(pow(float(obj1(0, j, k)) - float(obj2(0, j, k)), 2));

  float norm = sqrt(std::accumulate(coeff_diff.begin(), coeff_diff.end(), 0.0));

  return norm <= tolerance;

}

/**
 * This function defines the default test for mr_transform.
 * @param[in] Image file name.
 * @param[in] Expected result file name.
 * @param[in] Transform type.
 * @param[in] Number of scales.
 * @param[in] Filter type.
 * @return True or false.
 */
bool test_transform(char* image_file_name, char* result_file_name,
                   type_transform Transform, int Nbr_Plan,
                   type_undec_filter U_Filter){

  // Read in image file
  Ifloat image;
  io_read_ima_float(image_file_name, image);

  // Transform image and read expected result file
  MultiResol MR_Image, MR_Result;
  MR_Image.alloc (image.nl(), image.nc(), Nbr_Plan, Transform, NULL, NORM_L1,
                  -1, U_Filter);
  MR_Image.transform (image);
  MR_Result.read (result_file_name);

  // Test the difference between the transformation and the expectation
  return test_diff(MR_Image, MR_Result);

}

/**
 * This function defines the default test for mr_recons.
 * @param[in] Expected image file name.
 * @param[in] Transformed image file name.
 * @return True or false.
 */
bool test_recons(char* image_file_name, char* trans_file_name){

  // Read in image file
  Ifloat image;
  io_read_ima_float(image_file_name, image);

  // Read and reconstruct transformed image
  MultiResol MR_Image, MR_Trans;
  MR_Image.tabband() = &image;
  MR_Trans.read (trans_file_name);
  Ifloat rec (MR_Trans.size_ima_nl(), MR_Trans.size_ima_nc(), "Reconstruct.");
  MR_Trans.recons(rec);

  // Test the difference between the transformation and the expectation
  return test_diff(MR_Image, MR_Trans);

}

/*UNIT TESTS FOR MR_TRANSFORM*/
TEST_CASE("Unit tests on MR_TRANSFORM", "[mr_transform]") {

    REQUIRE(test_transform("../../tests/images/macos_cross.fits",
            "../../tests/images/macos_cross_t2_n4_U2.mr",
            TO_PAVE_BSPLINE, 4, DEF_UNDER_FILTER) == true);

    REQUIRE(test_transform("../../tests/images/macos_border.fits",
            "../../tests/images/macos_border_t16_n4_U2.mr",
            TO_PAVE_FEAUVEAU, 4, DEF_UNDER_FILTER) == true);

    REQUIRE(test_transform("../../tests/images/macos_circle.fits",
            "../../tests/images/macos_circle_t30_n4_U2.mr",
            TO_DIV_1, 4, DEF_UNDER_FILTER) == true);

}

/*UNIT TESTS FOR MR_RECONS*/
TEST_CASE("Unit tests on MR_RECONS", "[mr_recons]") {

    REQUIRE(test_recons("../../tests/images/macos_cross.fits",
            "../../tests/images/macos_cross_t2_n4_U2.mr") == true);

    REQUIRE(test_recons("../../tests/images/macos_border.fits",
            "../../tests/images/macos_border_t16_n4_U2.mr") == true);

    REQUIRE(test_recons("../../tests/images/macos_circle.fits",
            "../../tests/images/macos_circle_t30_n4_U2.mr") == true);

}
