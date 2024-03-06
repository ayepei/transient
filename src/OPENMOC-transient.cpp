#include "CPUSolver.h"
#include "log.h"
#include <array>
#include <iostream> 
#include "ThreadPrivateSolverTransient.h"
#include "Tcmfd.h"
#include "TransientSolver.h"
int main(int argc, char* argv[]) {
 /* Define simulation parameters */
  #ifdef OPENMP
  int num_threads = omp_get_num_procs();
  #else
  int num_threads = 1;
  #endif
  double track_spacing = 0.1;
  int num_azim = 4;
  double tolerance = 1e-5;
  int max_iters = 1000;
  double dt_cmfd = 1e-3;
  double dt_moc = 1e-2;
  /* Set logging information */
  setLogLevel("NORMAL");
  setOutputDirectory("/public/home/aicao/lsder/OpenMOC-transient/build/log");
  setLogfileName("2024.3.6");
  log_printf(TITLE, "Simulating the LRA Benchmark Problem...");

  /* Define material properties */
  log_printf(NORMAL, "Defining material properties...");

  const size_t num_groups = 2;
  std::map<std::string, std::array<double, num_groups> > nu_sigma_f;//Nu Fission XS
  std::map<std::string, std::array<double, num_groups> > sigma_f;//Fission XS
  std::map<std::string, std::array<double, num_groups*num_groups> > sigma_s;//Scattering XS
  std::map<std::string, std::array<double, num_groups> > chi;//CHI
  std::map<std::string, std::array<double, num_groups> > sigma_t;//Total XS
  std::map<std::string, std::array<double, num_groups> > sigma_a;     //Absorption XS
  std::map<std::string, std::array<double, num_groups> > Functional_Temperature;     //FunctionalTemperature = Absorption XS
  std::map<std::string, std::array<double, 1> > Temperature;     //Temperature
  std::map<std::string, bool> functional;     //Temperature
  std::map<std::string, std::array<double, 1> > FunctionalVariables;//  FunctionalVariables=Temperature
  std::map<std::string, std::array<double, num_groups> > DiffusionCoefficient;//瞬态独有(Ds)
  std::map<std::string, std::array<double, num_groups> > Buckling;//瞬态独有
  std::map<std::string, std::array<double, num_groups> > Gamma;//瞬态独有
  std::map<std::string, std::array<double, 3> > Time;//用于region_5
  std::map<std::string, std::array<double, num_groups> > FunctionalTime;     //FunctionalTime = Absorption XS 用于region_5
 


  /* Define region_1 cross-sections */
  nu_sigma_f["region_1"] = std::array<double, num_groups> {0.004602, 0.1091};
  sigma_f["region_1"] = std::array<double, num_groups> {0.002, 0.05};
  sigma_s["region_1"] = std::array<double, num_groups*num_groups>
      {0.232022, 0.02533, 0.00, 1.479479};
  chi["region_1"] = std::array<double, num_groups> {1.0, 0.0};
  sigma_t["region_1"] = std::array<double, num_groups> {0.2656, 1.5798};
  sigma_a["region_1"]=std::array<double, num_groups> {0.008252, 0.1003 };
  Functional_Temperature["region_1"]=sigma_a["region_1"];
  Temperature["region_1"]=std::array<double, 1> {300.0 };
  FunctionalVariables["region_1"]=Temperature["region_1"];
  DiffusionCoefficient["region_1"]=std::array<double,num_groups> {1.255, 0.211};
  Buckling["region_1"]=std::array<double,num_groups> {1e-4, 1e-4};
  Gamma["region_1"]=std::array<double,num_groups> {3.034e-3, 0.0};
  functional["region_1"]=false;

  /* Define region_2 cross-sections */
  //这里定义region_2-6 
  
  Temperature["region_2"]=std::array<double, 1> {300.0 };
  FunctionalVariables["region_2"]=Temperature["region_2"];
  sigma_a["region_2"]=std::array<double, num_groups> {0.007181, 0.07047 };
  Functional_Temperature["region_2"]=sigma_a["region_2"];
  sigma_t["region_2"] = std::array<double, num_groups> {0.2629, 1.7525};
  sigma_s["region_2"] = std::array<double, num_groups*num_groups>
      {0.22803, 0.02767, 0.00, 1.682071};
  sigma_f["region_2"] = std::array<double, num_groups> {0.002, 0.045};
  nu_sigma_f["region_2"] = std::array<double, num_groups> {0.004609, 0.08675};
  chi["region_2"] = std::array<double, num_groups> {1.0, 0.0};
  DiffusionCoefficient["region_2"]=std::array<double,num_groups> {1.268, 0.1902};
  Buckling["region_2"]=std::array<double,num_groups> {1e-4, 1e-4};
  Gamma["region_2"]=std::array<double,num_groups> {3.034e-3, 0.0};
  functional["region_2"]=false;

  /* Define region_3 cross-sections */
  Temperature["region_3"]=std::array<double, 1> {300.0 };
  FunctionalVariables["region_3"]=Temperature["region_3"];
  sigma_a["region_3"]=std::array<double, num_groups> {0.008002, 0.08344 };
  Functional_Temperature["region_3"]=sigma_a["region_3"];
  sigma_t["region_3"] = std::array<double, num_groups> {0.2648, 1.5941};
  sigma_s["region_3"] = std::array<double, num_groups*num_groups>
      {0.230588, 0.02617, 0.00, 1.510694};
  sigma_f["region_3"] = std::array<double, num_groups> {0.002, 0.045};
  nu_sigma_f["region_3"] = std::array<double, num_groups> {0.004663, 0.1021};
  chi["region_3"] = std::array<double, num_groups> {1.0, 0.0};
  DiffusionCoefficient["region_3"]=std::array<double,num_groups> {1.259, 0.2091};
  Buckling["region_3"]=std::array<double,num_groups> {1e-4, 1e-4};
  Gamma["region_3"]=std::array<double,num_groups> {3.034e-3, 0.0};
  functional["region_3"]=false;

  /* Define region_4 cross-sections */
  Temperature["region_4"]=std::array<double, 1> {300.0 };
  FunctionalVariables["region_4"]=Temperature["region_4"];
  sigma_a["region_4"]=std::array<double, num_groups> {0.008002, 0.073324 };
  Functional_Temperature["region_4"]=sigma_a["region_4"];
  sigma_t["region_4"] = std::array<double, num_groups> {0.2648, 1.5941};
  sigma_s["region_4"] = std::array<double, num_groups*num_groups>
      {0.230588, 0.02617, 0.00, 1.520810};
  sigma_f["region_4"] = std::array<double, num_groups> {0.002, 0.045};
  nu_sigma_f["region_4"] = std::array<double, num_groups> {0.004663, 0.1021};
  chi["region_4"] = std::array<double, num_groups> {1.0, 0.0};
  DiffusionCoefficient["region_4"]=std::array<double,num_groups> {1.259, 0.2091};
  Buckling["region_4"]=std::array<double,num_groups> {1e-4, 1e-4};
  Gamma["region_4"]=std::array<double,num_groups> {3.034e-3, 0.0};
  functional["region_4"]=false;

  /* Define region_5 cross-sections */
  //region_5为functional_material,在这里先正常赋值，在初始化的时候加入判断
  Temperature["region_5"]=std::array<double, 1> {300.0 };
  Time["region_5"]=std::array<double, 3> {0.0, 2.0, 3.0};
  FunctionalVariables["region_5"]=Temperature["region_5"];
  //  FunctionalVariables["region_5"]=Temperature["region_5"]+Time["region_5"];
//   FunctionalVariables_region5 = 

  sigma_a["region_5"]=std::array<double, num_groups> {0.008002, 0.08344 };
  // sigma_a["region_5"]=std::array<double, num_groups> {{0.008002, 0.08344}, {0.008002, 0.073324}, {0.008002, 0.073324}};
  double region5[3][2]={{0.008002, 0.08344}, {0.008002, 0.073324}, {0.008002, 0.073324}};
  
  Functional_Temperature["region_5"]=sigma_a["region_5"];
  FunctionalTime["region_5"]=sigma_a["region_5"];  
  sigma_t["region_5"] = std::array<double, num_groups> {0.2648, 1.5941};
  sigma_s["region_5"] = std::array<double, num_groups*num_groups>
      {0.230588, 0.02617, 0.00, 1.510694};
  sigma_f["region_5"] = std::array<double, num_groups> {0.002, 0.045};
  nu_sigma_f["region_5"] = std::array<double, num_groups> {0.004663, 0.1021};
  chi["region_5"] = std::array<double, num_groups> {1.0, 0.0};
  DiffusionCoefficient["region_5"]=std::array<double,num_groups> {1.259, 0.2091};
  Buckling["region_5"]=std::array<double,num_groups> {1e-4, 1e-4};
  Gamma["region_5"]=std::array<double,num_groups> {3.034e-3, 0.0};
  functional["region_5"]=true;

  /* Define region_6 cross-sections */
  sigma_a["region_6"]=std::array<double, num_groups> {0.0006034, 0.01911 };
  sigma_t["region_6"] = std::array<double, num_groups> {0.2652, 2.0938};
  sigma_s["region_6"] = std::array<double, num_groups*num_groups> 
      {0.217039, 0.04754, 0.00, 2.074692};
  sigma_f["region_6"] = std::array<double, num_groups> {0.0, 0.0};
  nu_sigma_f["region_6"] = std::array<double, num_groups> {0.0, 0.0};
  chi["region_6"] = std::array<double, num_groups> {1.0, 0.0};
  DiffusionCoefficient["region_6"]=std::array<double,num_groups> {1.257, 0.1592};
  Buckling["region_6"]=std::array<double,num_groups> {1e-4, 1e-4};
  functional["region_6"]=false;



  /* Create materials */
  log_printf(NORMAL, "Creating materials...");
  std::map<std::string, Material*> materials;

  std::map<std::string, std::array<double, num_groups> >::iterator it;
  int id_num = 0;
  /*遍历设置material*/
  for (it = sigma_t.begin(); it != sigma_t.end(); it++) {
    std::string name = it->first;


    if (Time.count(name) != 0){//find region5
        // std::cout<<name<<std::endl;

        materials[name] = reinterpret_cast<Material*>(new FunctionalMaterial(id_num));
        static_cast<FunctionalMaterial*>(materials[name])->setNumEnergyGroups(num_groups,Time[name].size());
        static_cast<FunctionalMaterial*>(materials[name])->setGamma(Gamma[name].data(),num_groups);
        static_cast<FunctionalMaterial*>(materials[name])->sigmaAFuncTemp(true);
        static_cast<FunctionalMaterial*>(materials[name])->setTime(Time[name].data(),Time[name].size());
        static_cast<FunctionalMaterial*>(materials[name])->sigmaAFuncTime(true);
        static_cast<FunctionalMaterial*>(materials[name])->setSigmaATime(Time[name].size(),num_groups,sigma_a[name].data());

    }
    else{
        if(Gamma.count(name)!=0){//find region1234
            // std::cout<<name<<std::endl;
            materials[name] = reinterpret_cast<Material*>(new FunctionalMaterial(id_num));
            static_cast<FunctionalMaterial*>(materials[name])->setNumEnergyGroups(num_groups,1);
            static_cast<FunctionalMaterial*>(materials[name])->setGamma(Gamma[name].data(),num_groups);
            static_cast<FunctionalMaterial*>(materials[name])->sigmaAFuncTemp(true);
            static_cast<FunctionalMaterial*>(materials[name])->setSigmaA(sigma_a[name].data(),num_groups);
        }

        else{//find region6
            // std::cout<<name<<std::endl;

            materials[name] = new Material(id_num);
            materials[name]->setNumEnergyGroups(num_groups);
            materials[name]->setSigmaA(sigma_a[name].data(),num_groups);
        }

    }

        materials[name]->setDifCoef(DiffusionCoefficient[name].data(),num_groups);//DiffusionCoefficient
        materials[name]->setBuckling(Buckling[name].data(),num_groups);//buckling
        materials[name]->setSigmaS(sigma_s[name].data(),num_groups*num_groups);//sigma_s
        materials[name]->setSigmaF(sigma_f[name].data(),num_groups);//sigma_f
        materials[name]->setNuSigmaF(nu_sigma_f[name].data(), num_groups);//nu_sigma_f
        materials[name]->setChi(chi[name].data(),num_groups);//chi
        id_num++;

        materials[name]->computeSigmaT();
        // std::cout<<"finish "<<name<<" set_value"<<std::endl;
    }


  /* Create surfaces */
    XPlane planes0(-82.5);
    XPlane planes1(82.5);
    YPlane planes2(-82.5);
    YPlane planes3(82.5);  

    planes0.setBoundaryType(REFLECTIVE);
    planes1.setBoundaryType(VACUUM);
    planes2.setBoundaryType(REFLECTIVE);
    planes3.setBoundaryType(VACUUM);  

    // std::cout<<"plane finished"<<std::endl;

//   /* Create circles for the fuel as well as to discretize the moderator into rings */

//     Circle fuel_radius(0.0, 0.0, 0.54);
//     Circle moderator_inner_radius(0.0, 0.0, 0.58);
//     Circle moderator_outer_radius(0.0, 0.0, 0.62);
//     std::cout<<"Circle finished"<<std::endl;

  /* Create cells and universes */
    log_printf(NORMAL, "Creating cells...");
    CellBasic* cell0 = new CellBasic(1,  materials["region_1"]->getId());
    CellBasic* cell1 = new CellBasic(2,  materials["region_2"]->getId());
    CellBasic* cell2 = new CellBasic(3,  materials["region_3"]->getId());
    CellBasic* cell3 = new CellBasic(4,  materials["region_4"]->getId());
    CellBasic* cell4 = new CellBasic(5,  materials["region_5"]->getId());
    CellBasic* cell5 = new CellBasic(6,  materials["region_6"]->getId());

    CellFill* cell6 = new CellFill(21, 31);
    CellFill* cell7 = new CellFill(22, 32);
    CellFill* cell8 = new CellFill(23, 33);
    CellFill* cell9 = new CellFill(24, 34);
    CellFill* cell10 = new CellFill(25, 35);
    CellFill* cell11 = new CellFill(26, 36);
    CellFill* cell12 = new CellFill(0, 7);

    cell12->addSurface(+1, &planes0);
    cell12->addSurface(-1,&planes1);
    cell12->addSurface(+1, &planes2);
    cell12->addSurface(-1,&planes3);

    log_printf(NORMAL, "Creating lattices...");

//define lattice assembly1
    Lattice* assembly1 = new Lattice(31,1.5,1.5);
    int mold1[10*10] =  {  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    assembly1->setLatticeCells(10,10,mold1);


 //define lattice assembly2
    Lattice* assembly2 = new Lattice(32,1.5,1.5);
    int mold2[10*10] =  {  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
    assembly2->setLatticeCells(10,10,mold2);   

 //define lattice assembly3
    Lattice* assembly3 = new Lattice(33,1.5,1.5);
    int mold3[10*10] =  {  3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3};
    assembly3->setLatticeCells(10,10,mold3);   

 //define lattice assembly4
    Lattice* assembly4 = new Lattice(34,1.5,1.5);
    int mold4[10*10] =  {  4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
    assembly4->setLatticeCells(10,10,mold4);   

 //define lattice assembly5
    Lattice* assembly5 = new Lattice(35,1.5,1.5);
    int mold5[10*10] =  {  5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
    assembly5->setLatticeCells(10,10,mold5);    
    
    
//define lattice assembly6
    Lattice* assembly6 = new Lattice(36,1.5,1.5);
    int mold6[10*10] =  {  6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
                           6, 6, 6, 6, 6, 6, 6, 6, 6, 6};
    assembly6->setLatticeCells(10,10,mold6);   


//define lattice core
    Lattice* core = new Lattice(7,15.0,15.0);
    int mold7[11*11] =  {26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
                         26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
                         23, 23, 23, 23, 23, 23, 23, 26, 26, 26, 26,
                         23, 23, 23, 23, 23, 23, 23, 24, 26, 26, 26,
                         22, 21, 21, 21, 21, 22, 22, 25, 25, 26, 26,
                         22, 21, 21, 21, 21, 22, 22, 25, 25, 26, 26,
                         21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26,
                         21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26,
                         21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26,
                         21, 21, 21, 21, 21, 21, 21, 23, 23, 26, 26,
                         22, 21, 21, 21, 21, 22, 22, 23, 23, 26, 26};
    core->setLatticeCells(11,11,mold7);   


//define mesh and geometry
    log_printf(NORMAL, "Creating cmfd mesh..");

    Mesh* mesh = new Mesh(MOC,true);


    log_printf(NORMAL, "Creating Geometry..");
    Geometry* geometry = new Geometry(mesh);

  for (it = sigma_t.begin(); it != sigma_t.end(); it++) {
        std::string name = it->first;
        geometry->addMaterial(materials[name]);
    }

    geometry->addCell(cell0);
    geometry->addCell(cell1);
    geometry->addCell(cell2);
    geometry->addCell(cell3);
    geometry->addCell(cell4);
    geometry->addCell(cell5);
    geometry->addCell(cell6);
    geometry->addCell(cell7);
    geometry->addCell(cell8);
    geometry->addCell(cell9);
    geometry->addCell(cell10);
    geometry->addCell(cell11);
    geometry->addCell(cell12);

    geometry->addLattice(assembly1);
    geometry->addLattice(assembly2);
    geometry->addLattice(assembly3);
    geometry->addLattice(assembly4);
    geometry->addLattice(assembly5);
    geometry->addLattice(assembly6);
    geometry->addLattice(core);

    geometry->initializeFlatSourceRegions();
    // int ii =1;
    // while(ii)
    //     printf("%d",ii);

//特征线生成
    TrackGenerator* track_generator = new TrackGenerator(geometry,num_azim,track_spacing);

    track_generator->generateTracks();

    log_printf(NORMAL, "Creating cmfd..");

    Cmfd* cmfd = new Cmfd(geometry);
    cmfd->setOmega(1.5);


//稳态求解器(瞬态版)
    ThreadPrivateSolverTransient* solver = new ThreadPrivateSolverTransient(geometry,track_generator,cmfd);
    solver->setNumThreads(num_threads);
    solver->setSourceConvergenceThreshold(tolerance);
    solver->initialize();//此处omp有问题，已经注释掉checkTrackSpacing中的omp部分


//Tcmfd

    Tcmfd* tcmfd = new Tcmfd(geometry);
    tcmfd->setOmega(1.5);
    double Matrix_lambda[2] ={0.0654, 1.35 };
    double Matrix_beta[2]={0.0054, 0.001087};    
    double Matrix_velocity[2]= {3e7, 3e5};
    tcmfd->setBeta(Matrix_beta,num_groups);
    tcmfd->setLambda(Matrix_lambda,num_groups);
    tcmfd->setVelocity(Matrix_velocity,num_groups);

    TransientSolver* transientSolver = new TransientSolver(geometry,tcmfd,cmfd,solver);
    transientSolver->setKappa(3.204e-11);
    transientSolver->setAlpha(3.83e-11);
    transientSolver->setNu(2.43);
    transientSolver->setDtMOC(dt_moc);
    transientSolver->setDtCMFD(dt_cmfd);
    transientSolver->setStartTime(0.0);
    transientSolver->setEndTime(3.0);
    transientSolver->setNumDelayGroups(num_groups);
    transientSolver->setTransientMethod("MAF");//有问题，但直接设值
    transientSolver->setPowerInit(1.e-6);

    transientSolver->solveInitialState();
    for(int i =0;i<int(1.5/dt_moc);i++)
        transientSolver->solveOuterStep();
    std::cout<<"chushihuawancheng"<<std::endl;
    log_printf(NORMAL, "Finished");



}
