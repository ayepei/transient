/*
 * TransientSolver.h
 *
 *  Created on: Dec 15, 2012
 *      Author: samuelshaner
 */

#ifndef TRANSIENTSOLVER_H_
#define TRANSIENTSOLVER_H_

#include <math.h>
#include <string>
#include <sstream>
#include <vector>
#include <ctime>
#include <time.h>
#include "Geometry.h"
#include "Quadrature.h"
#include "log.h"
#include "Tcmfd.h"
#include "Mesh.h"
#include "Cmfd.h"
#include "ThreadPrivateSolverTransient.h"
#include "Timer.h"
#include "TimeStepper.h"

class TransientSolver {
private:
    Geometry* _geom;
    Tcmfd* _tcmfd;
    Cmfd* _cmfd;
    Solver* _solver;
    TimeStepper* _ts;
    Mesh* _mesh;
    Mesh* _geom_mesh;
    Timer* _timer;

    double _dt_moc; // MOC方法的时间步长
    double _dt_cmfd;    // CMFD方法的时间步长
    double _power_init; // 初始功率水平
    double _vol_core; // 堆芯总体积
    double _k_eff_0; // 初始k有效值
    double _start_time;
    double _end_time;
    double _power_factor;            // 功率归一化因子
    solveType _solve_method;         // 求解方法类型(MOC/扩散)
    transientType _transient_method; // 瞬态方法类型(绝热/非绝热)
    int _ndg; // 缓发中子群数
    int _ng; // 能群数

    // double _nu;
    double _kappa; // 裂变能量释放系数(J/fission)
    double _alpha; // 温度反馈系数

    std::vector<double> _temp_core;   // 堆芯平均温度历史
    std::vector<double> _power_core;  // 堆芯功率历史
    std::vector<double> _time_tcmfd;  // TCMFD计算时间点
    std::vector<double> _temp_peak;   // 峰值温度历史
    std::vector<double> _power_peak;  // 峰值功率历史
    std::vector<double> _time_moc;    // MOC计算时间点
    std::vector<int> _num_amp_solves; // 振幅求解次数记录
    std::vector<double> _reactivity; // 反应性历史
    std::vector<int> _moc_iters; // MOC迭代次数历史
    std::string _log_file; // 日志文件路径

    double _temp_peak_value;  // 当前峰值温度
    double _power_peak_value; // 当前峰值功率

    /* global counters */
    int _amp_solve_counter; // 振幅求解计数器

    /* prolongation flag*/
    bool _prolongation; // 是否使用延拓标志

  public:
    TransientSolver(Geometry* geom, Tcmfd* tcmfd, Cmfd* cmfd,
		    Solver* solver=NULL);
    virtual ~TransientSolver();
    
    /* worker functions */
    double computeCoreTemp();  
    void computeVolCore();
    void sync(materialState state);
    void syncMaterials(materialState state);
    void copyPrecConc(materialState state_from, materialState state_to);
    void copyTemperature(materialState state_from, materialState state_to);
    void copyFieldVariables(materialState state_from, materialState state_to);
    double computePower(materialState state=CURRENT);
    void updateTemperatures();
    double computeResidual(solveType solve_type);
    void trimVectors(int len_conv);
    void initializePrecursorConc();
    // void updatePrecursorConc();
    void updatePrecursorConc(materialState state_from=PREVIOUS, materialState state_to=CURRENT);
    void initializeTimeStepper();
    void initializeTransientMaterials();
    void initializeTransientLogfile();
    void solveInitialState();
    void solveOuterStep();
    void mapPrecConc();
    void logStep();
    
    /* setters */
    void setDtMOC(double dt);
    void setDtCMFD(double dt);  
    void setKappa(double kappa);
    // void setNu(double nu);
    void setAlpha(double alpha);
    void setNumDelayGroups(int num_groups);
    void setTransientMethod(const char* trans_type);
    void setPowerInit(double power);
    void setStartTime(double time);
    void setEndTime(double time);
    double getPower();
    double getTime();
    double getTemp();
    void setProlongation(bool prolong);



};

#endif /* TRANSIENTSOLVER_H_ */
