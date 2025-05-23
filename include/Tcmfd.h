/**
 * @file Tcmfd.h
 * @brief The Tcmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef TCMFD_H_
#define TCMFD_H_

#ifdef __cplusplus
#define _USE_MATH_DEFINES
#include <math.h>
#include "Material.h"
#include "TimeStepper.h"
// #include "FunctionalMaterial.h"
#include "Geometry.h"
#include "Timer.h"
#include "Mesh.h"
#include "linalg_functions.h"
#include <vector>
#endif


class Tcmfd{
protected:
  TimeStepper* _time_stepper;
  Mesh* _mesh;
  Geometry* _geom;
  Timer* _timer;

  /* integer values */
  int _cells_x; // x方向网格数
  int _cells_y;
  int _num_groups;
  int _nc; // 总单元数(cells_x * cells_y * _num_groups)
  int _num_delay_groups; // 缓发中子组数

  /* float values */
  double _omega; // 松弛因子
  double _k_eff_0; // 初始k有效值
  double _beta_sum; // 缓发中子份额总和
  double _dt_cmfd; // CMFD时间步长
  double _conv_criteria; // 收敛判据
  
  /* matrix and vector objects */
  double* _b_prime; // 修正的右端项
  double* _b; // 右端项
  double* _A; // 损失矩阵
  double* _M; // 裂变矩阵
  double* _phi_old; // 旧通量
  double* _phi_new;
  double* _phi_temp;
  double* _snew; // 新源项
  double* _sold;
  double* _AM; // A和M的乘积
 
  /* materials parameters */
  double* _beta; // 缓发中子份额
  double* _lambda; // 缓发中子衰变常数
  double* _velocity; // 中子速度

  /* solve method (DIFFUSION or MOC) */
  solveType _solve_method; 
  transientType _transient_method;
  bool _initial_state; 
 
  
public:
  Tcmfd(Geometry* geometry, double criteria=1e-8);
  virtual ~Tcmfd();
  void constructMatrices(bool frequency=false);
  void solveTCMFD();
  void setTransientType(transientType trans_type);
  void setKeff0(double keff_0);
  transientType getTransientType();
  void setTimeStepper(TimeStepper* _time_stepper);
  
  void setBeta(double* beta, int num_delay_groups);
  void setLambda(double* decay_const, int num_delay_groups);
  void setVelocity(double* velocity, int num_groups);
  void setDtCMFD(double dt);

  double* getBeta();
  double* getLambda();
  double* getVelocity();
  
  void setInitialState(bool state);
  void setNumDelayGroups(int num_groups);
  int getNumDelayGroups();

  Mesh* getMesh();
  void setOmega(double omega);

  void checkNeutronBalance();
  void checkNeutronBalance2();
  void computeFrequency();
};

#endif /* TCMFD_H_ */
