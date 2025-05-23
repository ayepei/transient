/**
 * @file Cmfd.h
 * @brief The Cmfd class.
 * @date October 14, 2013
 * @author Sam Shaner, MIT, Course 22 (shaner@mit.edu)
 */

#ifndef MESH_H_
#define MESH_H_

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <utility>
#include <sstream>
#include "log.h"
#include "LocalCoords.h"
#include "Surface.h"
#include "Material.h"
#include "Quadrature.h"
#include "TimeStepper.h"
#include "Timer.h"
#include "linalg_functions.h"

/**
 * Solve types
 */
enum solveType {
	DIFFUSION,
	MOC
};

enum transientType {
  NONE,
  THETA,
  MAF,
  ADIABATIC
};


class Mesh {

private:

  /* physical mesh size */
  double _length_x; // 网格在x方向的总长度
  double _length_y;

  /* pointer to quadrature object */
  Quadrature* _quad; // 积分对象
  TimeStepper* _ts;  // 时间步进器
  Timer* _timer;     // 计时器

  /* cmfd level CMFD网格级别*/
  int _mesh_level;

  /* number of cells in x and y directions */
  int _cells_x; // x方向的网格单元数量
  int _cells_y;

  /* number of groups */
  int _num_groups;        // 能量群数量
  int _num_delay_groups;  // 缓发中子群数量

  /* number of surface current values 表面流数量*/
  int _num_currents;

  /* number of fsrs 精细空间区域(FSR)数量*/
  int _num_fsrs;

  /* number of azim angles 方位角数量*/
  int _num_azim;    

  /* array of boundary enums 边界类型*/
  boundaryType* _boundaries;

  /* array of mesh cell volumes 网格单元体积*/
  double* _volumes;

  /* array of mesh surface currents 表面流数据*/
  double* _currents;

  /* vector of vectors of fsrs in each mesh cell 每个网格单元包含的FSR*/
  std::vector< std::vector<int> > _cell_fsrs;

  /* cmfd and acceleration flag */
  bool _cmfd_on;
  bool _acceleration;

  /* relaxation factor on d_tilde 松弛因子*/
  double _relax_factor;  

  /* map of fluxes mapped onto mesh 不同材料状态下的通量*/
  std::map<materialState, double*> _fluxes;
  // std::map<materialState, double*> _currents;
  // std::map<materialState, double*> _frequency;
  
  /* frequency array 存储频率相关的数据*/
  double* _frequency;

  /* map of fsrs to cells FSR到网格单元的映射*/
  int* _FSRs_to_cells;

  /* materials array 材料数组*/
  Material** _materials;

  /* array of fsr bounds 存储FSR的索引信息*/
  int* _fsr_indices;

  /* array of lenghts of each mesh cell in x and y directions */
  double* _lengths_x;  // 每个网格单元在x方向的长度
  double* _lengths_y;

  /* array of cell bounds in x and y direction */
  double* _bounds_x; // x方向的网格边界坐标
  double* _bounds_y;

  Material** _FSR_materials;  // FSR对应的材料
  FP_PRECISION* _FSR_volumes; // 存储每个精细空间区域(FSR)的体积
  FP_PRECISION* _FSR_fluxes; // FSR通量

  /* bool to toggle optically thick diffusion correction factor 是否使用光学厚扩散修正*/
  bool _optically_thick;

  /* solve method (DIFFUSION or MOC) */
  solveType _solve_method;
  transientType _transient_method;

  double* _lambda;     // 衰变常数
  double* _beta;       // 缓发中子份额
  double _beta_sum;    // 总缓发中子份额
  double* _velocity;   // 中子速度
  bool _initial_state; // 是否为初始状态
  double _k_eff_0;     // 初始有效增殖因子
  double _dt_moc;      // MOC时间步长
  
public:
  Mesh(solveType solve_type=MOC, bool cmfd_on=false,
       double relax_factor=0.6, int mesh_level=-1);
  virtual ~Mesh();
  void initialize();
  void setFSRBounds();
  void setCellBounds();

  // void computeDs(double relax_factor=-1.0, materialState state=FSR_OLD);
  // void computeXS(Mesh* mesh=NULL, materialState state=FSR);
  void computeDs(double relax_factor=-1.0, materialState state=SHAPE);
  void computeXS(Mesh* mesh=NULL, materialState state=SHAPE);
  double computeDiffCorrect(double d, double h);
  void updateMOCFlux();

  /* get mesh parameters */
  double getLengthX();
  double getLengthY();
  int getCellsX();
  int getCellsY();
  int getNumCells();
  boundaryType getBoundary(int side);
  int getNumCurrents();
  double getFlux(int cell_id, int group, materialState state=CURRENT);
  std::vector<std::vector<int> >* getCellFSRs();
  Material** getMaterials();
  double* getVolumes();
  double* getFluxes(materialState state);
  // double* getFrequencies(materialState state);
  double* getFrequency();
  double* getLengthsX();
  double* getLengthsY();
  // double* getCurrents(materialState state);
  double* getCurrents();
  int getMeshLevel();
  
  /* set mesh parameters */
  void setLengthX(double length_x);
  void setLengthY(double length_y);
  void setCellLengthX(int cell_num, double length_x);
  void setCellLengthY(int cell_num, double length_y);
  void setCellsX(int cells_x);
  void setCellsY(int cells_y);
  void setSurfaceCurrents(double* surface_currents);
  void setVolume(double volume, int cell_num);
  void setMeshLevel(int cmfd_level);
  void setFlux(materialState state, double* flux);
  void eraseFlux(materialState state);

  /* set general problem specs */
  void setNumGroups(int num_groups);
  void setNumAzim(int num_azim);
  void setNumFSRs(int num_fsrs);
  void setAcceleration(bool accel);
  void setOpticallyThick(bool thick);
  void setRelaxFactor(double relax_factor);

  /* get generation problem specs */
  int getNumGroups();
  int getNumFSRs();
  bool getCmfdOn();
  bool getAcceleration();
  bool getOpticallyThick();
  double getRelaxFactor();
  solveType getSolveType();

  /* worker functions */
  int findMeshCell(double x, double y);
  int findMeshSurface(int fsr_id, LocalCoords* coord, int angle);
  void printCurrents();
  void splitCorners();
  void setBoundary(int side, boundaryType boundary);
  int getCellNext(int cell_num, int surface_id);
  int findCellId(LocalCoords* coord);
  void initializeMaterials(std::map<int, Material*>* materials, int* fsrs_to_mats);
  void initializeSurfaceCurrents();
  void createNewFlux(materialState state);
  // void createNewFrequency(materialState state);
  // void createNewCurrent(materialState state);
  void copyFlux(materialState from_state, materialState to_state);
  // void copyFrequency(materialState from_state, materialState to_state);
  // void copyCurrent(materialState from_state, materialState to_state);
  void copyDs(materialState from_state, materialState to_state);
  void dumpFlux(materialState state);
  void dumpXS();

  void setFSRMaterials(Material** FSR_materials);
  void setFSRVolumes(FP_PRECISION* FSR_volumes);
  void setFSRFluxes(FP_PRECISION* scalar_flux);
  Material** getFSRMaterials();
  FP_PRECISION* getFSRFluxes();

  void geomSetMaterials(Material** FSR_materials);
  void geomSetVolumes(FP_PRECISION* FSR_volumes);

  void setNumDelayGroups(int num_groups);
  int getNumDelayGroups();

  void setLambda(double* decay_constant, int ndg);
  void setBeta(double* beta, int ndg);
  void setVelocity(double* velocity, int ng);
  double* getLambda();
  double* getBeta();
  double getBetaSum();
  double* getVelocity();
  transientType getTransientType();
  void setTransientType(transientType trans_method);
  void setInitialState(bool init);
  bool getInitialState();
  void setKeff0(double k_eff_0);
  double getKeff0();
  void setDtMOC(double dt);
  double getDtMOC();
  void setFSRToCell(int fsr_id, int cell_id);
  int* getFSRsToCells();

  void reconstructFineFlux(double* geom_shape, double* mesh_flux);
  void computeFineShape(double* geom_shape, double* mesh_flux);
  // void interpolateCurrent(double ratio);
  void interpolateDs(double ratio);
  void interpolateFlux(double ratio);
  void normalizeFlux(double scale_val);
  void normalizeDs(double scale_val);
  void zeroDs();

  void setTimeStepper(TimeStepper* ts);
  TimeStepper* getTimeStepper();
  double getPrecursorConc(int fsr_id, int dg);
  double getSigmaA(int fsr_id, int g);
  void dumpFSRs();
  double getTemperature(int fsr_id);
};

#endif /* MESH_H_ */
