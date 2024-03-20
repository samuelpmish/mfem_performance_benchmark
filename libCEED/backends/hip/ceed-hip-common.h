// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and other CEED contributors.
// All Rights Reserved. See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-2-Clause
//
// This file is part of CEED:  http://github.com/ceed

#ifndef CEED_COMMON_HIP_H
#define CEED_COMMON_HIP_H

#include <ceed.h>
#include <ceed/backend.h>
#include <hip/hip_runtime.h>
#if (HIP_VERSION >= 50200000)
#include <hipblas/hipblas.h>  // IWYU pragma: export
#else
#include <hipblas.h>  // IWYU pragma: export
#endif

#define QUOTE(...) #__VA_ARGS__

#define CeedChk_Hip(ceed, x)                             \
  do {                                                   \
    hipError_t hip_result = x;                           \
    if (hip_result != hipSuccess) {                      \
      const char *msg = hipGetErrorName(hip_result);     \
      return CeedError((ceed), CEED_ERROR_BACKEND, msg); \
    }                                                    \
  } while (0)

#define CeedChk_Hipblas(ceed, x)                             \
  do {                                                       \
    hipblasStatus_t hipblas_result = x;                      \
    if (hipblas_result != HIPBLAS_STATUS_SUCCESS) {          \
      const char *msg = hipblasGetErrorName(hipblas_result); \
      return CeedError((ceed), CEED_ERROR_BACKEND, msg);     \
    }                                                        \
  } while (0)

#define CeedCallHip(ceed, ...)        \
  do {                                \
    hipError_t ierr_q_ = __VA_ARGS__; \
    CeedChk_Hip(ceed, ierr_q_);       \
  } while (0)

#define CeedCallHipblas(ceed, ...)         \
  do {                                     \
    hipblasStatus_t ierr_q_ = __VA_ARGS__; \
    CeedChk_Hipblas(ceed, ierr_q_);        \
  } while (0)

#define CASE(name) \
  case name:       \
    return #name
// LCOV_EXCL_START
CEED_UNUSED static const char *hipblasGetErrorName(hipblasStatus_t error) {
  switch (error) {
    CASE(HIPBLAS_STATUS_SUCCESS);
    CASE(HIPBLAS_STATUS_NOT_INITIALIZED);
    CASE(HIPBLAS_STATUS_ALLOC_FAILED);
    CASE(HIPBLAS_STATUS_INVALID_VALUE);
    CASE(HIPBLAS_STATUS_ARCH_MISMATCH);
    CASE(HIPBLAS_STATUS_MAPPING_ERROR);
    CASE(HIPBLAS_STATUS_EXECUTION_FAILED);
    CASE(HIPBLAS_STATUS_INTERNAL_ERROR);
    default:
      return "HIPBLAS_STATUS_UNKNOWN_ERROR";
  }
}
// LCOV_EXCL_STOP

typedef struct {
  int                    device_id;
  hipblasHandle_t        hipblas_handle;
  struct hipDeviceProp_t device_prop;
  int                    opt_block_size;
} Ceed_Hip;

CEED_INTERN int CeedInit_Hip(Ceed ceed, const char *resource);

CEED_INTERN int CeedDestroy_Hip(Ceed ceed);

#endif  // CEED_COMMON_HIP_H
