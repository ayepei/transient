#include "Material.h"

int Material::_n = 0;


static int auto_id = 10000;


/**
 * @brief Returns an auto-generated unique material ID.
 * @details This method is intended as a utility method for user's writing
 *          OpenMOC input files. The method makes use of a static material
 *          ID which is incremented each time the method is called to enable
 *          unique generation of monotonically increasing IDs. The method's
 *          first ID begins at 10000. Hence, user-defined material IDs greater
 *          than or equal to 10000 is prohibited.
 */
int material_id() {
    int id = auto_id;
    auto_id++;
    return id;
}


/**
 * @brief Constructor sets the ID and unique ID for the material.
 * @param id the user-defined id for the material
 */
Material::Material(short int id) {

    _uid = _n;
    _id = id;
    _n++;
    _type = BASE;

    _sigma_t = NULL;
    _sigma_a = NULL;
    _sigma_a_ref = NULL;
    _sigma_s = NULL;
    _sigma_s_ref = NULL;
    _sigma_f = NULL;
    _nu_sigma_f = NULL;
    _chi = NULL;
    _dif_coef = NULL;
    _dif_hat = NULL;
    _dif_tilde = NULL;
    _buckling = NULL;

    _fissionable = false;
    _data_aligned = false;
    
    _temperature = new double[6];
    
    for (int i = 0; i < 6; i++)
      _temperature[i] = 300.0;

    return;
}


/**
 * @brief Destructor deletes all cross-section data.
 */
Material::~Material() { 

    /* If data is vector aligned */
    if (_data_aligned) {
        if (_sigma_t != NULL)
	    _mm_free(_sigma_t);

	if (_sigma_a != NULL)
	    _mm_free(_sigma_a);

	if (_sigma_s != NULL)
	    _mm_free(_sigma_s);
	
	if (_sigma_f != NULL)
	    _mm_free(_sigma_f);
	
	if (_nu_sigma_f != NULL)
	    _mm_free(_nu_sigma_f);
	
	if (_chi != NULL)
	    _mm_free(_chi);

	if (_dif_coef != NULL)
	    delete [] _dif_coef;

	if (_dif_hat != NULL)
	    delete [] _dif_hat;

	if (_dif_tilde != NULL)
	    delete [] _dif_tilde;

	if (_buckling != NULL)
	    delete [] _buckling;

    }

    /* Data is not vector aligned */
    else {
        if (_sigma_t != NULL)
	    delete [] _sigma_t;

	if (_sigma_a != NULL)
	    delete [] _sigma_a;

	if (_sigma_s != NULL)
	    delete [] _sigma_s;
	
	if (_sigma_f != NULL)
	    delete [] _sigma_f;
	
	if (_nu_sigma_f != NULL)
	    delete [] _nu_sigma_f;
	
	if (_chi != NULL)
	    delete [] _chi;

	if (_dif_coef != NULL)
	    delete [] _dif_coef;

	if (_dif_hat != NULL)
	    delete [] _dif_hat;

	if (_dif_tilde != NULL)
	    delete [] _dif_tilde;

	if (_buckling != NULL)
	    delete [] _buckling;
    }
}


/**
 * @brief Return the material's unique ID.
 * @return the material's unique ID
 */
int Material::getUid() const {
    return _uid;
}

/**
 * @brief Return the material's user-defined ID
 * @return the material's user-defined ID
 */
short int Material::getId() const {
    return _id;
}


/**
 * @brief Returns the number of energy groups for this material's nuclear data.
 * @return the number of energy groups
 */
int Material::getNumEnergyGroups() const {
    return _num_groups;
}


/**
 * @brief Return the array of the material's total cross-sections.
 * @return the pointer to the material's array of total cross-sections
 */
double* Material::getSigmaT() {
    if (_sigma_t == NULL)
        log_printf(ERROR, "Unable to return material %d's total "
                        "cross-section since it has not yet been set", _id);
    return _sigma_t;
}



/**
 * @brief Return the array of the material's absorption cross-sections.
 * @return the pointer to the material's array of absorption cross-sections
 */
double* Material::getSigmaA() {
    if (_sigma_a == NULL)
        log_printf(ERROR, "Unable to return material %d's absorption "
                        "cross-section since it has not yet been set", _id);

    return _sigma_a;
}


/**
 * @brief Return the array of the material's scattering cross-section matrix.
 * @return the pointer to the material's array of scattering cross-sections
 */
double* Material::getSigmaS() {
    if (_sigma_s == NULL)
        log_printf(ERROR, "Unable to return material %d's scattering "
                        "cross-section since it has not yet been set", _id);

    return _sigma_s;
}

/**
 * @brief Return the array of the material's reference scattering cross-section matrix.
 * @return the pointer to the material's array of reference scattering cross-sections
 */
double* Material::getSigmaSRef() {
    if (_sigma_s_ref == NULL)
        log_printf(ERROR, "Unable to return material %d's scattering "
                        "cross-section since it has not yet been set", _id);

    return _sigma_s_ref;
}

/**
 * @brief Return the array of the material's fission cross-sections.
 * @return the pointer to the material's array of fission cross-sections
 */
double* Material::getSigmaF() {

    if (_sigma_f == NULL)
        log_printf(ERROR, "Unable to return material %d's fission "
                        "cross-section since it has not yet been set", _id);

    return _sigma_f;
}


/**
 * @brief Return the array of the material's fission cross-sections
 *        multiplied by nu \f$ \nu \f$.
 * @return the pointer to the material's array of fission cross-sections 
 *         multiplied by nu \f$ \nu \f$
 */
double* Material::getNuSigmaF() {
    if (_nu_sigma_f == NULL)
        log_printf(ERROR, "Unable to return material %d's nu times fission "
                        "cross-section since it has not yet been set", _id);

    return _nu_sigma_f;
}


/**
 * @brief Return the array of the material's chi \f$ \chi \f$.
 * @return the pointer to the material's array of chi \f$ \chi \f$ values
 */
double* Material::getChi() {
    if (_chi == NULL)
        log_printf(ERROR, "Unable to return material %d's chi spectrum "
                        "since it has not yet been set", _id);

    return _chi;
}


/**
 * @brief Return the array of the material's diffusion coefficients.
 * @return the pointer to the material's array of diffusion coefficients
 */
double* Material::getDifCoef() {

  if (_dif_coef == NULL){
    _dif_coef = new double[_num_groups];
    for (int e = 0; e < _num_groups; e++)
      _dif_coef[e] = 0.0;
  }

    return _dif_coef;
}


/**
 * @brief Return the array of the material's surface diffusion coefficients.
 * @return the pointer to the material's array of surface diffusion coefficients
 */
double* Material::getDifHat() {

  if (_dif_hat == NULL){
    _dif_hat = new double[4*_num_groups];
    for (int e = 0; e < 4*_num_groups; e++)
      _dif_hat[e] = 0.0;
  }


    return _dif_hat;
}

/**
 * @brief Return the array of the material's CMFD correction to the
 * surface diffusion coefficients.
 * @return the pointer to the material's array of CMFD correction to
 * the surface diffusion coefficients
 */
double* Material::getDifTilde() {

  if (_dif_tilde == NULL){
    _dif_tilde = new double[6*4*_num_groups];
    for (int e = 0; e < 6*4*_num_groups; e++)
    // _dif_tilde = new double[4*_num_groups];
    // for (int e = 0; e < 4*_num_groups; e++)
      _dif_tilde[e] = 0.0;
  }

    return _dif_tilde;
}


/**
 * @brief Return the array of the material's CMFD correction to the
 * surface diffusion coefficients.
 * @return the pointer to the material's array of CMFD correction to
 * the surface diffusion coefficients
 */
double* Material::getBuckling() {

  if (_buckling == NULL){
    _buckling = new double[_num_groups];
    for (int e = 0; e < _num_groups; e++)
      _buckling[e] = 0.0;
  }

    return _buckling;
}


/**
 * @brief Returns whether or not the material contains a fissionable (non-zero)
 *        fission cross-section.
 * @return true if fissionable, false otherwise
 */
bool Material::isFissionable() {
    return _fissionable;
}


/**
 * @brief Returns true if the data is vector aligned, false otherwise (default).
 * @return Whether or not the materials data is vector aligned
 */
bool Material::isDataAligned() {
    return _data_aligned;
}


/**
 * @brief Returns the rounded up number of energy groups to fill an integral
 *        number of vector lengths
 * @return The number of vector-aligned energy groups
 */
int Material::getNumVectorGroups() {
    return _num_vector_groups;
}



/**
 * @brief Set the number of energy groups for this material.
 * @param num_groups the number of energy groups.
 */
void Material::setNumEnergyGroups(const int num_groups) {

    if (num_groups < 0)
        log_printf(ERROR, "Unable to set the number of energy groups for "
                   "material %d to %d", _num_groups);

    _num_groups = num_groups;

    /* Free old memory arrays if they have already been allocated 
     * for a previous simulation */
    /* If data is vector aligned */
    if (_data_aligned) {
        if (_sigma_t != NULL)
	    _mm_free(_sigma_t);

	if (_sigma_a != NULL)
	    _mm_free(_sigma_a);

	if (_sigma_s != NULL)
	    _mm_free(_sigma_s);
	
	if (_sigma_f != NULL)
	    _mm_free(_sigma_f);
	
	if (_nu_sigma_f != NULL)
	    _mm_free(_nu_sigma_f);
	
	if (_chi != NULL)
	    _mm_free(_chi);
    }

    /* Data is not vector aligned */
    else {
        if (_sigma_t != NULL)
	    delete [] _sigma_t;

	if (_sigma_a != NULL)
	    delete [] _sigma_a;

	if (_sigma_s != NULL)
	    delete [] _sigma_s;
	
	if (_sigma_f != NULL)
	    delete [] _sigma_f;
	
	if (_nu_sigma_f != NULL)
	    delete [] _nu_sigma_f;
	
	if (_chi != NULL)
	    delete [] _chi;

	if (_dif_coef != NULL)
	    delete [] _dif_coef;

	if (_dif_hat != NULL)
	    delete [] _dif_hat;

	if (_dif_tilde != NULL)
	    delete [] _dif_tilde;

	if (_buckling != NULL)
	    delete [] _buckling;
    }

    _sigma_t = new double[_num_groups];
    _sigma_a = new double[_num_groups];
    _sigma_f = new double[_num_groups];
    _nu_sigma_f = new double[_num_groups];
    _chi = new double[_num_groups];
    _sigma_s = new double[_num_groups*_num_groups];
    _buckling = new double[_num_groups];

    for (int g = 0; g < _num_groups; g++)
      _buckling[g] = 0.0;
}


/**
 * @brief Set the material's array of total cross-sections.
 * @param xs the array of total cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaT(double* xs, int num_groups) {
    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set sigma_t with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++)
        _sigma_t[i] = xs[i];
}


void Material::computeSigmaT(){

    /* set sigma t to sigma_a */
    for (int i=0; i < _num_groups; i++)
        _sigma_t[i] = _sigma_a[i];
    
    /* add scattering xs */
    for (int e=0; e < _num_groups; e++){
	for (int g=0; g < _num_groups; g++)
	    _sigma_t[e] += _sigma_s[g*_num_groups+e];
    }
}

/**
 * @brief Set the material's total cross-section for some energy group.
 * @param xs the total cross-section (\f$ \Sigma_t [cm^1] \f$)
 * @param group the energy group
 */
void Material::setSigmaTByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_t for group %d for material "
                 "%d which contains %d energy groups", group, _uid, _num_groups);

    _sigma_t[group] = xs;
}


/**
 * @brief Set the material's array of absorption cross-sections.
 * @details This method is intended to be called from 
 * @param xs the array of absorption cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaA(double* xs, int num_groups) {
    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set sigma_a with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i = 0; i < _num_groups; i++)
      _sigma_a[i] = xs[i];

    if (_dif_coef != NULL && _buckling != NULL){
      for (int i = 0; i < _num_groups; i++)
	_sigma_a[i] += _dif_coef[i] * _buckling[i];
    }
}


/**
 * @brief Set the material's absorption cross-section for some energy group.
 * @param xs the absorption cross-section (\f$ \Sigma_a [cm^1] \f$)
 * @param group the energy group
 */
void Material::setSigmaAByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_a for group %d for material "
                 "%d which contains %d energy groups", group, _uid, _num_groups);

    _sigma_a[group] = xs;
}




/**
 * @brief Set the material's 2D array of scattering cross-sections. 
 * @details This assumes that the scattering matrix passed in has the standard 
 *          notation: the ij element is for scattering from group i to j. For 
 *          efficient caching of the elements of this matrix during fixed 
 *          source iteration, the matrix transpose is what is actually stored 
 *          in the material
 * @param xs the array of scattering cross-sections
 * @param num_groups_squared the number of energy groups squared
 */
void Material::setSigmaS(double* xs, int num_groups_squared) {
 
    if (_num_groups*_num_groups != num_groups_squared)
        log_printf(ERROR, "Unable to set sigma_s with %f groups for material "
		   "%d which contains %d energy groups", 
		  //  float(sqrt(num_groups_squared)), _num_groups);
      float(sqrt(num_groups_squared)), _id, _num_groups);


    for (int i=0; i < _num_groups; i++) {
      for (int j=0; j < _num_groups; j++){
	_sigma_s[j*_num_groups+i] = xs[i*_num_groups+j];
      }   
    }
}


/**
 * @brief Set the material's scattering cross-section for some energy group.
 * @param xs the scattering cross-section (\f$ \Sigma_s [cm^1] \f$)
 * @param group1 the row index in the scattering matrix
 * @param group2 the column index in the scattering matrix
 */
void Material::setSigmaSByGroup(double xs, int group1, int group2) {

    if (group1 < 0 || group2 < 0 || group1 >= _num_groups || group2 >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_s for group %d,%d for material "
                 "%d which contains %d energy groups", group1, group2, _uid,
                 _num_groups);

    _sigma_s[_num_groups*group1 + group2] = xs;
}



/**
 * @brief Set the material's array of fission cross-sections.
 * @param xs the array of fission cross-sections
 * @param num_groups the number of energy groups
 */
void Material::setSigmaF(double* xs, int num_groups) {
    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++)
        _sigma_f[i] = xs[i];

    /* Determine whether or not this material is fissionable */
    _fissionable = false;

    for (int i=0; i < _num_groups; i++) {
        if (_sigma_f[i] > 0.0) {
	    _fissionable = true;
	    return;
	}
    }
}


/**
 * @brief Set the material's fission cross-section for some energy group.
 * @param xs the fission cross-section (\f$ \Sigma_f [cm^1] \f$)
 * @param group the energy group
 */
void Material::setSigmaFByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set sigma_f for group %d for material "
                 "%d which contains %d energy groups", group, _uid, _num_groups);

    _sigma_f[group] = xs;

    /* Determine whether or not this material is fissionable */
    _fissionable = false;

    for (int i=0; i < _num_groups; i++) {
        if (_sigma_f[i] > 0.0) {
	    _fissionable = true;
	    return;
	}
    }
}


/**
 * @brief Set the material's array of fission cross-sections multiplied by
 *         \f$ \nu \f$
 * @param xs the array of fission cross-sections multiplied by nu 
 *        \f$ \nu \f$
 * @param num_groups the number of energy groups 
*/
void Material::setNuSigmaF(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set nu_sigma_f with %d groups for material "
                 "%d which contains %d energy groups", num_groups, _uid, _num_groups);

    for (int i=0; i < _num_groups; i++)
        _nu_sigma_f[i] = xs[i];
}


/**
 * @brief Set the material's fission cross-section multiplied by \f$ \nu \f$
 *        for some energy group.
 * @param xs the fission cross-section (\f$ \nu\Sigma_f [cm^1] \f$)
 * @param group the energy group
 */
void Material::setNuSigmaFByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
      log_printf(ERROR, "Unable to set nu_sigma_f for group %d for material "
                 "%d which contains %d energy groups", group, _uid);

    _nu_sigma_f[group] = xs;
}



/**
 * @brief Set the material's array of \f$ \chi \f$ values.
 * @param xs the array of chi \f$ \chi \f$ values
 * @param num_groups the number of energy groups 
 */
void Material::setChi(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set chi with %d groups for material "
                 "%d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++)
        _chi[i] = xs[i];
}


/**
 * @brief Set the material's chi value for some energy group.
 * @param xs the chi value (\f$ \Chi \f$)
 * @param group the energy group
 */
void Material::setChiByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set chi for group %d for material "
		   "%d which contains %d energy groups", 
		   group, _num_groups, _uid);

    _chi[group] = xs;
}


/**
 * @brief Set the material's array of diffusion coefficients.
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setDifCoef(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
    		  "for material %d which contains %d energy groups", num_groups,
                 _num_groups);

    if (_dif_coef == NULL)
      _dif_coef = new double[_num_groups];

    for (int i=0; i < _num_groups; i++)
        _dif_coef[i] = xs[i];
}



/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDifCoefByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    if (_dif_coef == NULL){
      _dif_coef = new double[_num_groups];
      for (int i=0; i < _num_groups; i++)
        _dif_coef[i] = 0.0;
    }

    _dif_coef[group] = xs;
}


/**
 * @brief Set the material's array of diffusion coefficients.
 * @param xs the array of diffusion coefficents
 * @param num_groups the number of energy groups
 */
void Material::setBuckling(double* xs, int num_groups) {

    if (_num_groups != num_groups)
      log_printf(ERROR, "Unable to set diffusion coefficient with %d groups "
    		  "for material %d which contains %d energy groups", num_groups,
                 _num_groups);

    for (int i=0; i < _num_groups; i++)
        _buckling[i] = xs[i];
}


/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setBucklingByGroup(double xs, int group) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    _buckling[group] = xs;
}


/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
void Material::setDifHatByGroup(double xs, int group, int surface) {

    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    if (_dif_hat == NULL){
      _dif_hat = new double[4*_num_groups];
      for (int i=0; i < _num_groups*4; i++)
        _dif_hat[i] = 0.0;
    }

    _dif_hat[surface*_num_groups + group] = xs;
}



/**
 * @brief Set the material's diffusion coefficient for some energy group.
 * @param xs the diffusion coefficient
 * @param group the energy group
 */
// void Material::setDifTildeByGroup(double xs, int group, int surface) {
  void Material::setDifTildeByGroup(double xs, int group, int surface, materialState state) {


    if (group < 0 || group >= _num_groups)
        log_printf(ERROR, "Unable to set diffusion coefficient correction for group %d"
        		" for material %d which contains %d energy groups",
		   group, _num_groups, _uid);

    int ngs = 4*_num_groups;

    if (_dif_tilde == NULL){
      _dif_tilde = new double[6*ngs];
      for (int i=0; i < 6*ngs; i++)
      // _dif_tilde = new double[ngs];
      // for (int i=0; i < ngs; i++)
	      _dif_tilde[i] = 0.0;
    }
    _dif_tilde[int(state)*ngs + surface*_num_groups + group] = xs;
    // _dif_tilde[surface*_num_groups + group] = xs;
}


/**
 * @brief Checks if the total cross-section for this material is equal to the
 *        absorption plus scattering cross-sections for all energy groups.
 * @details If the total cross-section does not equal the absorption plus 
 *          scattering cross-section within SIGMA_T_THRESH (defined in
 *          openmoc/src/host/configurations.h) then this method exits OpenMOC.
 */
void Material::checkSigmaT() {

    if (_num_groups == 0)
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since the number of energy groups has not been set", _id);
    if (_sigma_t == NULL) 
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since its total cross-section has not been set", _id);
    if (_sigma_a == NULL) 
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since its absorption cross-section has not been set", _id);
    if (_sigma_s == NULL)
        log_printf(ERROR, "Unable to verify material %d's total cross-section "
                "since its scattering cross-section has not been set", _id);

    double calc_sigma_t;

    /* Loop over all energy groups */
    for (int i=0; i < _num_groups; i++) {

        /* Initialize the calculated total xs to the absorption xs */
        if (_type == BASE)
	  calc_sigma_t = _sigma_a[i];
	else
	  calc_sigma_t = _sigma_a_ref[i];

        /* Increment calculated total xs by scatter xs for each energy group */
        // for (int j=0; j < _num_groups; j++)
          for (int j=0; j < _num_groups; j++){
	          if (_type == BASE){
            calc_sigma_t += _sigma_s[i+j*_num_groups];
	      }
	  else{
            calc_sigma_t += _sigma_s_ref[i+j*_num_groups];
	  } 
	}
        /* Check if the calculated and total match up to certain threshold */
        if (fabs(calc_sigma_t - _sigma_t[i]) > SIGMA_T_THRESH) {
            log_printf(ERROR, "Material id = %d has a different total "
                       "cross-section than the sum of its scattering and "
                       "absorption cross-sections for group %d: "
                       "sigma_t = %f, calc_sigma_t = %f", _id, i, _sigma_t[i],
                       calc_sigma_t);
        }
    }
    
    return;
}


/**
 * @brief Converts this material's attributes to a character array 
 *        representation.
 * @details The character array returned includes the user-defined ID,
 *          and each of the absorption, total, fission, nu multiplied by
 *          fission and scattering cross-sections and chi for all energy
 *          groups.
 * @return character array of this member's attributes
 */
std::string Material::toString() {

    std::stringstream string;

    string << "Material id = " << _id;

    if (_sigma_a != NULL) {
        string << "\n\t\tSigma_a = ";
        for (int e = 0; e < _num_groups; e++)
            string << _sigma_a[e] << ", ";
    }

    if (_sigma_t != NULL) {
        string << "\n\t\tSigma_t = ";
        for (int e = 0; e < _num_groups; e++)
            string << _sigma_t[e] << ", ";
    }

    if (_sigma_f != NULL) {
        string << "\n\t\tSigma_f = ";
        for (int e = 0; e < _num_groups; e++)
            string << _sigma_f[e] << ", ";
    }


    if (_nu_sigma_f != NULL) {
        string << "\n\t\tnu_sigma_f = ";
        for (int e = 0; e < _num_groups; e++)
            string << _nu_sigma_f[e] << ", ";
    }

    if (_sigma_s != NULL) {
        string << "\n\t\tSigma_s = \n\t\t";
        for (int G = 0; G < _num_groups; G++) {
  	    for (int g = 0; g < _num_groups; g++)
                string << _sigma_s[G+g*_num_groups] << "\t\t ";
            string << "\n\t\t";
        }
    }

    if (_chi != NULL) {
        string << "Chi = ";
        for (int e = 0; e < _num_groups; e++)
            string << _chi[e] << ", ";
    }

    if (_dif_coef != NULL) {
        string << "Diffusion Coefficient = ";
        for (int e = 0; e < _num_groups; e++)
            string << _dif_coef[e] << ", ";
    }

    if (_buckling != NULL) {
        string << "Buckling = ";
        for (int e = 0; e < _num_groups; e++)
            string << _buckling[e] << ", ";
    }

    return string.str();
}


/**
 * @brief Prints a string representation of all of the material's objects to
 *        the console.
 */
void Material::printString() {
    printf("%s", toString().c_str());
}


/**
 * @brief Aligns the cross-section data structures 
 */
void Material::alignData() {

    if (_data_aligned)
        return;

    if (_num_groups <= 0)
        log_printf(ERROR, "Unable to align material %d data since the "
		   "cross-sections have not yet been set\n", _uid);

    _num_vector_groups = (_num_groups / VEC_LENGTH) + 1;

    /* Allocate memory for the new aligned xs data */
    int size = _num_vector_groups * VEC_LENGTH * sizeof(double);

    double* new_sigma_t = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_sigma_a = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_sigma_f = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_nu_sigma_f = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    double* new_chi = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    
    /* The scattering matrix will be the number of vector groups 
     * wide (SIMD) and the actual number of groups long since 
     * instructions are not SIMD in this dimension */
    size = _num_vector_groups * VEC_LENGTH * _num_vector_groups * VEC_LENGTH * sizeof(double);
    double* new_sigma_s = (double*)_mm_malloc(size, VEC_ALIGNMENT);
    
    /* Initialize data structures to ones for sigma_t since it is used to
     * divide the source in the solver, and zeroes for everything else */
    size = _num_vector_groups * VEC_LENGTH * sizeof(double);
    for (int i=0; i < _num_vector_groups * VEC_LENGTH; i++) {
        new_sigma_t[i] = 1.0;
	new_sigma_a[i] = 0.0;
	new_sigma_f[i] = 0.0;
	new_nu_sigma_f[i] = 0.0;
	new_chi[i] = 0.0;
    }
    
    size *= _num_vector_groups * VEC_LENGTH;
    memset(new_sigma_s, 0.0, size);
    
    /* Copy materials data from unaligned arrays into new aligned arrays */
    size = _num_groups * sizeof(double);
    memcpy(new_sigma_t, _sigma_t, size);
    memcpy(new_sigma_a, _sigma_a, size);
    memcpy(new_sigma_f, _sigma_f, size);
    memcpy(new_nu_sigma_f, _nu_sigma_f, size);
    memcpy(new_chi, _chi, size);

    for (int e=0; e < _num_groups; e++) {
        memcpy(new_sigma_s, _sigma_s, size);
	new_sigma_s += _num_vector_groups * VEC_LENGTH;
	_sigma_s += _num_groups;
    }

    _sigma_s -= _num_groups * _num_groups;
    
    /* Reset the new scattering cross section array pointer */
    new_sigma_s -= _num_vector_groups * VEC_LENGTH * _num_groups;

    /* Delete the old unaligned arrays */
    delete [] _sigma_t;
    delete [] _sigma_a;
    delete [] _sigma_f;
    delete [] _nu_sigma_f;
    delete [] _chi;
    delete [] _sigma_s;
    	
    /* Set the material's array pointers to the new aligned arrays */
    _sigma_t = new_sigma_t;
    _sigma_a = new_sigma_a;
    _sigma_f = new_sigma_f;
    _nu_sigma_f = new_nu_sigma_f;
    _chi = new_chi;
    _sigma_s = new_sigma_s;
    
    _data_aligned = true;

    return;
}


Material* Material::clone(){

  Material* material_clone = new Material(getId());
  
  material_clone->setNumEnergyGroups(_num_groups);
  material_clone->setSigmaT(_sigma_t, _num_groups);
  material_clone->setSigmaF(_sigma_f, _num_groups);
  material_clone->setNuSigmaF(_nu_sigma_f, _num_groups);
  material_clone->setChi(_chi, _num_groups);
  
  if (_dif_coef != NULL)
    material_clone->setDifCoef(_dif_coef, _num_groups);  
  
  if (_buckling != NULL)
    material_clone->setBuckling(_buckling, _num_groups);  

  if (_dif_hat != NULL)
    for (int i = 0; i < _num_groups; i++)
      for (int s = 0; s < 4; s++)
	material_clone->setDifHatByGroup(_dif_hat[s*_num_groups+i], i, s);  

  int ngs = _num_groups*4;

  if (_dif_tilde != NULL){
      for (int i = 0; i < _num_groups; i++){
	  for (int s = 0; s < 4; s++){
	      material_clone->setDifTildeByGroup(_dif_tilde[0*ngs + s*_num_groups+i], i, s, PREVIOUS_CONV);  
	      material_clone->setDifTildeByGroup(_dif_tilde[1*ngs + s*_num_groups+i], i, s, PREVIOUS);  
	      material_clone->setDifTildeByGroup(_dif_tilde[2*ngs + s*_num_groups+i], i, s, CURRENT);  
	      material_clone->setDifTildeByGroup(_dif_tilde[3*ngs + s*_num_groups+i], i, s, FORWARD);  
	      material_clone->setDifTildeByGroup(_dif_tilde[4*ngs + s*_num_groups+i], i, s, FORWARD_PREV);  
	      material_clone->setDifTildeByGroup(_dif_tilde[5*ngs + s*_num_groups+i], i, s, SHAPE);  

	      // material_clone->setDifTildeByGroup(_dif_tilde[s*_num_groups+i], i, s);  
	  }
      }
  }


  copySigmaA(material_clone);
  copySigmaS(material_clone);

  material_clone->setTemperature(PREVIOUS, getTemperature(PREVIOUS));
  material_clone->setTemperature(PREVIOUS_CONV, getTemperature(PREVIOUS_CONV));
  material_clone->setTemperature(CURRENT, getTemperature(CURRENT));
  material_clone->setTemperature(FORWARD, getTemperature(FORWARD));
  // material_clone->setTemperature(FSR, getTemperature(FSR));
  material_clone->setTemperature(FORWARD_PREV, getTemperature(FORWARD_PREV));
  material_clone->setTemperature(SHAPE, getTemperature(SHAPE));
  return material_clone;

}


void Material::copySigmaS(Material* material){

  for (int i=0; i < _num_groups; i++){
    for (int j=0; j < _num_groups; j++){
      material->setSigmaSByGroup(_sigma_s[i*_num_groups+j], i, j);
    }
  }
}


void Material::copySigmaA(Material* material){

  for (int i=0; i < _num_groups; i++){
      material->setSigmaAByGroup(_sigma_a[i], i);
  }
}

void Material::initializeTemperature(double temperature){

  _temperature[0] = temperature;
  _temperature[1] = temperature;
  _temperature[2] = temperature;
  _temperature[3] = temperature;
  _temperature[4] = temperature;
  _temperature[5] = temperature;
}


void Material::setTemperature(materialState state, double temperature){
  _temperature[(int)state] = temperature;
}


double Material::getTemperature(materialState state){
  return _temperature[(int)state];
}

void Material::copyTemperature(materialState state_from, materialState state_to){
  _temperature[(int)state_to] = _temperature[(int)state_from];
}

materialType Material::getType(){
  return _type;
}

/**
 * @brief Constructor sets the ID and unique ID for the material.
 * @param id the user-defined id for the material
 */
FunctionalMaterial::FunctionalMaterial(short int id) : Material(id){

  _time = NULL;
  _type = FUNCTIONAL;
  _ts = NULL;
  _gamma = NULL;
  
  _sigma_a_func_temp = false;  
  _sigma_a_func_time = false;
  _sigma_s_func_time = false;
  _conserve_sigma_t  = true;
  _sigma_a_ref = NULL;
  _sigma_s_ref = NULL;

  _temperature = new double[6];
  
  for (int i = 0; i < 6; i++)
    _temperature[i] = 300.0;
  
  return;
}


/**
 * @brief Destructor deletes all cross-section data.
 */
FunctionalMaterial::~FunctionalMaterial() {

  if (_sigma_a_ref != NULL)
    delete [] _sigma_a_ref;  
  if (_sigma_s_ref != NULL)
    delete [] _sigma_s_ref;    
}




/**
 * @brief Set the number of energy groups for this material.
 * @param num_groups the number of energy groups.
 */
void FunctionalMaterial::setNumEnergyGroups(const int num_groups, const int num_time_steps) {

  Material::setNumEnergyGroups(num_groups);
  
  log_printf(DEBUG, "Setting material's num energy groups");
  
  if (num_groups < 0)
    log_printf(ERROR, "Unable to set the number of energy groups for "
	       "material %d to %d", _num_groups);
  
  _num_time_steps = num_time_steps;

  _sigma_a_ref = new double[num_groups*num_time_steps];
  _sigma_s_ref = new double[num_groups*num_groups*num_time_steps];

  for (int i = 0; i < num_groups*num_time_steps; i++)
    _sigma_a_ref[i] = 0.0;

  for (int i = 0; i < num_groups*num_groups*num_time_steps; i++)
    _sigma_s_ref[i] = 0.0;
  _gamma = new double[num_groups];
}


/**
 * @brief Set the material's array of absorption scattering cross-sections.
 * @brief Set the material's array of absorption cross-sections.
 * @details This method is intended to be called from 
 * @param xs the array of absorption scattering cross-sections
 * @param xs the array of absorption cross-sections
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setSigmaA(double* xs, int num_groups) {

    if (_num_groups != num_groups)
	log_printf(ERROR, "Unable to set sigma_a with %d groups for material "
		   "%d which contains %d energy groups", num_groups,
		   _num_groups);
  
    Material::setSigmaA(xs, num_groups);

    for (int i=0; i < _num_groups; i++)
	_sigma_a_ref[i] = xs[i];
}


/** @brief Set the material's array of scattering cross-sections.
 * @details This method is intended to be called from 
 * @param xs the array of scattering cross-sections
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setSigmaS(double* xs, int num_groups_squared) {

    if (_num_groups*_num_groups != num_groups_squared)
	log_printf(ERROR, "Unable to set sigma_s with %d groups for material "
		   "%d which contains %d energy groups", float(sqrt(num_groups_squared)), _id,
		   _num_groups);

    Material::setSigmaS(xs, num_groups_squared);
    // 这里保存的时候使用的是转置，为什么
    for (int i=0; i < _num_groups; i++) {
      for (int j=0; j < _num_groups; j++){
	    _sigma_s_ref[j*_num_groups+i] = xs[i*_num_groups+j];
      }   
    }
}


/**
 * @brief Set the material's array of absorption scattering cross-sections.
 * @details This method is intended to be called from
 * @param xs the array of absorption scattering cross-sections
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setSigmaATime(int num_time_steps, int num_groups, double* xs) {

  if (_num_groups != num_groups)
    log_printf(ERROR, "Unable to set sigma_a with %d groups for material "
	       "%d which contains %d energy groups", num_groups,
	       _num_groups);

  /* load _sigma_t_ref with all xs */
  for (int i = 0; i < num_time_steps*num_groups; i++)
      _sigma_a_ref[i] = xs[i];
  
  for (int i = 0; i < _num_groups; i++)
      _sigma_a[i] = xs[i];

  if (_dif_coef != NULL && _buckling != NULL){
    for (int i = 0; i < _num_groups; i++)
      _sigma_a[i] += _dif_coef[i] * _buckling[i];
  }
}

/**
 * @brief Set the material's array of scattering cross-sections.
 * @details This method is intended to be called from
 * @param xs the array of scattering cross-sections
 * @param num_groups the number of energy groups
 */
void FunctionalMaterial::setSigmaSTime(int num_time_steps, int num_groups_squared, double* xs) {

  if (_num_groups != int(sqrt(num_groups_squared)))
    log_printf(ERROR, "Unable to set sigma_s with %d groups for material "
	       "%d which contains %d energy groups", int(sqrt(num_groups_squared)),
	       _num_groups);

  int ng = _num_groups;

  /* set the reference scattering xs array */
  for (int i = 0; i < num_time_steps; i++){
    for (int k=0; k < ng; k++) {
      for (int j=0; j < ng; j++){
	_sigma_s_ref[i*ng*ng + j*ng + k] = xs[i*ng*ng + k*ng + j];
      }   
    }
  }

  /* set the callable scattering xs array */
  for (int k=0; k < ng; k++) {
    for (int j=0; j < ng; j++){
      _sigma_s[j*ng + k] = xs[k*ng + j];
    }   
  }
}

FunctionalMaterial* FunctionalMaterial::clone(){

  FunctionalMaterial* to_mat = new FunctionalMaterial(getId());

  /* copy flags */
  to_mat->sigmaAFuncTime(_sigma_a_func_time);
  to_mat->sigmaAFuncTemp(_sigma_a_func_temp);
  to_mat->sigmaSFuncTime(_sigma_s_func_time);
  to_mat->setConserveSigmaT(_conserve_sigma_t);
  /* set num energy groups */
  to_mat->setNumEnergyGroups(_num_groups, _num_time_steps);

  /* copy xs */
  to_mat->setSigmaT(_sigma_t, _num_groups);
  copySigmaS(to_mat);
  copySigmaSRef(to_mat);
  to_mat->setSigmaF(_sigma_f, _num_groups);
  to_mat->setNuSigmaF(_nu_sigma_f, _num_groups);
  to_mat->setChi(_chi, _num_groups);
  
  if (_buckling != NULL)
    to_mat->setBuckling(_buckling, _num_groups);
				
  if (_dif_coef != NULL)
    to_mat->setDifCoef(_dif_coef, _num_groups);

  if (_sigma_a_func_time)
    to_mat->setSigmaATime(_num_time_steps, _num_groups, _sigma_a_ref);
  else
    to_mat->setSigmaA(_sigma_a_ref, _num_groups);
  
  to_mat->setTemperature(PREVIOUS, getTemperature(PREVIOUS));
  to_mat->setTemperature(PREVIOUS_CONV, getTemperature(PREVIOUS_CONV));
  to_mat->setTemperature(CURRENT, getTemperature(CURRENT));
  to_mat->setTemperature(FORWARD, getTemperature(FORWARD));
  // to_mat->setTemperature(FSR, getTemperature(FSR));
  to_mat->setTemperature(FORWARD_PREV, getTemperature(FORWARD_PREV));
  to_mat->setTemperature(SHAPE, getTemperature(SHAPE));

  if (_gamma != NULL)
    to_mat->setGamma(_gamma, _num_groups);
  
  if (_time != NULL)
    to_mat->setTime(_time, _num_time_steps);

  if (_ts != NULL)
    to_mat->setTimeStepper(_ts);

  return to_mat;
}


void FunctionalMaterial::setTime(double* time, int num_time_steps) {
  
  _time = new double[num_time_steps];
  
  for (int i=0; i < num_time_steps; i++)
    _time[i] = time[i];
}


double* FunctionalMaterial::getTime(){
  return _time;
}


void FunctionalMaterial::sigmaAFuncTemp(bool func_temp){
  _sigma_a_func_temp = func_temp;
}


void FunctionalMaterial::sigmaAFuncTime(bool func_time){
  _sigma_a_func_time = func_time;
}

void FunctionalMaterial::sigmaSFuncTime(bool func_time){
  _sigma_s_func_time = func_time;
}
// 这个函数的作用？确保了材料截面数据在不同时间步和温度状态下的正确更新，同时保持了物理守恒性，主要用于同步更新材料在不同状态下的截面数据
/* sync current cross sections with current time and temperature */
void FunctionalMaterial::sync(materialState state){
  
    double sigma_s_out;
    double sigma_s_group;

    /* SIGMA_A */
    for (int g = 0; g < _num_groups; g++){
    
        sigma_s_out = 0.0;
      
	if (_sigma_a_func_time)
	    _sigma_a[g] = interpolateXS(_sigma_a_ref, state, g);
	else
	    _sigma_a[g] = _sigma_a_ref[g];
	
	if (_sigma_a_func_temp)
	    _sigma_a[g] = _sigma_a[g] * (1.0 + _gamma[g] * 
	                  (pow(getTemperature(state),0.5) - pow(300.0, 0.5)));

	/* adjust self scattering to conserve total xs */
	// for (int G = 0; G < _num_groups; G++){
  if (_conserve_sigma_t){//为true

	  /* adjust self scattering to conserve total xs */
	  for (int G = 0; G < _num_groups; G++){
        if (G != g)
	        sigma_s_out += _sigma_s[G*_num_groups + g];
	}
      
	// _sigma_s[g*_num_groups + g] = 1.0 / (3.0 * _dif_coef[g]) - _sigma_a[g] - sigma_s_out;    

	// _sigma_a[g] += _dif_coef[g] * _buckling[g];

	// _sigma_t[g] = _sigma_a[g] + sigma_s_out + _sigma_s[g*_num_groups + g];
    _sigma_s[g*_num_groups + g] = 1.0 / (3.0 * _dif_coef[g]) - _sigma_a[g] - sigma_s_out;    

	  _sigma_a[g] += _dif_coef[g] * _buckling[g];

	  _sigma_t[g] = _sigma_a[g] + sigma_s_out + _sigma_s[g*_num_groups + g];
	}
	else{

	  sigma_s_group = 0.0;

	  /* sync scattering cross section */
	  if (_sigma_s_func_time){
	    for (int G = 0; G < _num_groups; G++){
	      _sigma_s[G*_num_groups + g] = interpolateScatterXS(_sigma_s_ref, state, g, G);
	      sigma_s_group += _sigma_s[G*_num_groups + g];
	    }
	  }
	  else{
	    for (int G = 0; G < _num_groups; G++){
	      sigma_s_group += _sigma_s[G*_num_groups + g];
	    }
	  }

	  /* recompute sigma_t */
	  _sigma_t[g] = _sigma_a[g] + sigma_s_group;
	}
    }
}


void FunctionalMaterial::setGamma(double* gamma, int num_groups){

  for (int g = 0; g < num_groups; g++)
    _gamma[g] = gamma[g];
}


double* FunctionalMaterial::getGamma(){
  return _gamma;
}

// 使用线性插值计算中间时刻的截面值
double FunctionalMaterial::interpolateXS(double* xs_ref, materialState state, int group){

  double time = _ts->getTime(state);
  double dt, dxs;
  double xs = xs_ref[group];
  
  for (int i = 1; i < _num_time_steps; i++){
    if (time < _time[i] + 1e-8){
      dt = _time[i] - _time[i - 1];
      // 计算截面变化量
      dxs = xs_ref[i*_num_groups + group] - xs_ref[(i-1)*_num_groups + group];
      xs = xs_ref[(i-1)*_num_groups + group] + (time - _time[i-1]) / dt * dxs;
      break;
    }
  }
  
  return xs;
}

// 散射矩阵插值
double FunctionalMaterial::interpolateScatterXS(double* xs_ref, materialState state, int group_from, int group_to){

  double time = _ts->getTime(state);
  double dt, dxs;
  int ng = _num_groups;
  double xs = xs_ref[group_to*ng + group_from];


  for (int i = 1; i < _num_time_steps; i++){
    if (time < _time[i] + 1e-8){
      dt = _time[i] - _time[i-1];
      dxs = xs_ref[i*ng*ng + group_to*ng + group_from] - xs_ref[(i-1)*ng*ng + group_to*ng + group_from];
      xs = xs_ref[(i-1)*ng*ng + group_to*ng + group_from] + (time - _time[i-1]) / dt * dxs;
      break;
    }
  }

  return xs;
}

void FunctionalMaterial::initializeTransientProps(double num_delay_groups, bool cmfd_mesh){

  _num_delay_groups = num_delay_groups;
  _prec_conc = new double[_num_delay_groups*6];
  _prec_freq = new double[_num_delay_groups*6];
  // _prec_conc = new double[_num_delay_groups*7];
  // _prec_freq = new double[_num_delay_groups*7];

}

void FunctionalMaterial::setPrecConc(materialState state, double conc, int group){
  _prec_conc[(int)state * _num_delay_groups + group] = conc;
}

void FunctionalMaterial::setPrecFreq(materialState state, double freq, int group){
  _prec_freq[(int)state *_num_delay_groups + group] = freq;
}

double FunctionalMaterial::getPrecConc(materialState state, int group){
  return _prec_conc[(int)state *_num_delay_groups + group];
}

double FunctionalMaterial::getPrecFreq(materialState state, int group){
  return _prec_freq[(int)state*_num_delay_groups + group];
}

void FunctionalMaterial::copyPrecConc(materialState state_from, materialState state_to){
  
  for (int dg = 0; dg < _num_delay_groups; dg++)
    _prec_conc[(int)state_to *_num_delay_groups + dg] = _prec_conc[(int)state_from *_num_delay_groups + dg];
}

void FunctionalMaterial::copyPrecFreq(materialState state_from, materialState state_to){
 
 for (int dg = 0; dg < _num_delay_groups; dg++)
    _prec_freq[(int)state_to *_num_delay_groups + dg] = _prec_freq[(int)state_from *_num_delay_groups + dg];
}


void FunctionalMaterial::setTimeStepper(TimeStepper* ts){
  _ts = ts;
}


// /**
//  * @brief Return the array of the material's absorption cross-sections.
//  * @return the pointer to the material's array of absorption cross-sections
//  */
// double FunctionalMaterial::getSigmaAByValue(materialState state, int group) {
void FunctionalMaterial::setConserveSigmaT(bool conserve_sigma_t){
  _conserve_sigma_t = conserve_sigma_t;
}

  // double xs = _sigma_a_ref[group];
  
  // if (_sigma_a_func_time)
  //   xs = interpolateXS(_sigma_a_ref, state, group);
  
  // if (_sigma_a_func_temp)
  //     xs = xs * (1.0 + _gamma[group] * (pow(getTemperature(state),0.5) - pow(300.0, 0.5)));
  
  // return xs;
  void FunctionalMaterial::copySigmaSRef(Material* material){

  double* xs_ref = material->getSigmaSRef();

  for (int i = 0; i < _num_time_steps*_num_groups*_num_groups; i++)
    xs_ref[i] = _sigma_s_ref[i];
}
