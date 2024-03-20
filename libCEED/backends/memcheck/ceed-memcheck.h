// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and other CEED contributors.
// All Rights Reserved. See the top-level LICENSE and NOTICE files for details.
//
// SPDX-License-Identifier: BSD-2-Clause
//
// This file is part of CEED:  http://github.com/ceed

#ifndef CEED_MEMCHECK_H
#define CEED_MEMCHECK_H

#include <ceed.h>
#include <ceed/backend.h>

typedef struct {
  int         mem_block_id;
  bool        is_write_only_access;
  CeedScalar *array;
  CeedScalar *array_allocated;
  CeedScalar *array_owned;
  CeedScalar *array_borrowed;
  CeedScalar *array_read_only_copy;
} CeedVector_Memcheck;

typedef struct {
  const CeedScalar **inputs;
  CeedScalar       **outputs;
  bool               setup_done;
} CeedQFunction_Memcheck;

typedef struct {
  int   mem_block_id;
  void *data;
  void *data_allocated;
  void *data_owned;
  void *data_borrowed;
  void *data_read_only_copy;
} CeedQFunctionContext_Memcheck;

CEED_INTERN int CeedVectorCreate_Memcheck(CeedSize n, CeedVector vec);

CEED_INTERN int CeedQFunctionCreate_Memcheck(CeedQFunction qf);

CEED_INTERN int CeedQFunctionContextCreate_Memcheck(CeedQFunctionContext ctx);

#endif  // CEED_MEMCHECK_H
