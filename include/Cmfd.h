/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef CMFD_H_
#define CMFD_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <utility>
#include <math.h>
#include <unordered_map>
#include <limits.h>
#include <string>
#include <sstream>
#include <queue>
#include <iostream>
#include <fstream>
#include "Quadrature.h"
#include "log.h"
#include "Mesh.h"
#include "Material.h"
#include "Surface.h"
#include "Geometry.h"
#include "Timer.h"
#include "linalg_functions.h"
#include <omp.h>
#endif

class Cmfd {
protected:

  Geometry* _geometry;
  Mesh* _mesh;
  Timer* _timer;  

  /* keff */
  double _k_eff;

  /* matrix and vector objects */
  double* _A; // 损失矩阵
  double* _M; // 裂变矩阵
  double* _phi_old; // 旧通量
  double* _phi_new; // 新通量
  double* _sold; // 旧源项
  double* _snew; // 新源项
  double* _phi_temp; // 临时通量
  double* _b; // 右端项
  double *_b_prime; // 修正的右端项 


  /* float values */
  // double _conv_criteria;
  double _conv_linear; // 线性迭代收敛准则
  double _conv_nonlinear; // 非线性迭代收敛准则
  double _omega; // 松弛因子

  /* integer values */
  int _cx; // x方向网格数
  int _cy; // y方向网格数
  int _ng; // 能群数
  int _nc; // 总单元数(cx*cy*_ng)

  /* solve method (DIFFUSION or MOC) */
  solveType _solve_method;
  
  double* _AM; // A和M的乘积
  double* _y;  // 中间计算结果


public:
	
  // Cmfd(Geometry* geometry, double criteria=1e-8);
  Cmfd(Geometry* geometry, double conv_linear=1.e-8, double conv_nonlinear=1.e-6);

  virtual ~Cmfd();

  virtual void constructMatrices();
  double computeKeff();
  void rescaleFlux();
  double* getA();
  double* getM();
  Mesh* getMesh();
  double getKeff();
  void setOmega(double omega);
  void checkNeutronBalance();
  void setNumThreads(int num_threads);
};

#endif /* CMFD_H_ */
