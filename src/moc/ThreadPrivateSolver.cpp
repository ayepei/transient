#include "ThreadPrivateSolver.h"


/**
 * @brief Constructor initializes array pointers for tracks and materials.
 * @details The constructor retrieves the number of energy groups and flat
 *          source regions and azimuthal angles from the geometry and track
 *          generator. The constructor initalizes the number of threads to a 
 *          default of 1.
 * @param geometry an optional pointer to the geometry
 * @param track_generator an optional pointer to the trackgenerator
 */
ThreadPrivateSolver::ThreadPrivateSolver(Geometry* geometry, 
					 TrackGenerator* track_generator,
					 Cmfd* cmfd) :
  CPUSolver(geometry, track_generator, cmfd) {

    _thread_flux = NULL;
    _thread_currents = NULL;
}


/**
 * @brief Destructor calls Solver subclass destructor to deletes arrays
 *        for fluxes and sources.
 */
ThreadPrivateSolver::~ThreadPrivateSolver() { 

    if (_thread_flux != NULL) {
        delete [] _thread_flux;
	_thread_flux = NULL;
    }

    if (_thread_currents != NULL) {
        delete [] _thread_currents;
	_thread_currents = NULL;
    }

}


/**
 * @brief Allocates memory for track boundary angular fluxes and leakages
 *        flat source region scalar fluxes.
 * @details Deletes memory for old flux arrays if they were allocated from
 *          previous simulation.
 */
void ThreadPrivateSolver::initializeFluxArrays() {

    CPUSolver::initializeFluxArrays();
   
    /* Delete old flux arrays if they exist */
    if (_thread_flux != NULL)
        delete [] _thread_flux;

    /* Delete old current arrays if they exist */
    if (_thread_currents != NULL)
        delete [] _thread_currents;

    int size;

    /* Allocate memory for the flux and leakage arrays */
    try{
	/* Allocate a thread local array of FSR scalar fluxes */
	size = _num_threads * _num_FSRs * _num_groups * sizeof(FP_PRECISION);
	_thread_flux = new FP_PRECISION[size];

	/* Allocate a thread local array of mesh cell surface currents */
	if (_cmfd->getMesh()->getCmfdOn()){ 
	  size = _num_threads * _num_mesh_cells * 8 * _num_groups * sizeof(double);
	  _thread_currents = new double[size];
	}
    }
    catch(std::exception &e) {
        log_printf(ERROR, "Could not allocate memory for the solver's fluxes. "
		   "Backtrace:%s", e.what());
    }
}


/**
 * @brief Set the scalar flux for each energy group inside each flat source 
 *        region to a constant value.
 * @details This method also flattens the thread private flat source region
 *          scalar flux array.
 * @param value the value to assign to each flat source region flux
 */
void ThreadPrivateSolver::flattenFSRFluxes(FP_PRECISION value) {

    CPUSolver::flattenFSRFluxes(value);

    /* Flatten the thread private flat source region scalar flux array */
    #pragma omp parallel for schedule(guided)
    for (int tid=0; tid < _num_threads; tid++) {
        for (int r=0; r < _num_FSRs; r++) {
	    for (int e=0; e < _num_groups; e++) {
	        _thread_flux(tid,r,e) = 0.0;
	    }
        }
    }

    return;
}


/**
 * @brief Set the surface currents for each energy group inside each
 * 			mesh cell to zero
 */
void ThreadPrivateSolver::zeroSurfaceCurrents() {

	CPUSolver::zeroSurfaceCurrents();

    #pragma omp parallel for schedule(guided)
	for (int tid=0; tid < _num_threads; tid++){
		for (int r=0; r < _num_mesh_cells; r++) {
			for (int s=0; s < 8; s++) {
				for (int e=0; e < _num_groups; e++)
					_thread_currents(tid,r*8+s,e) = 0.0;
			}
		}
	}

    return;
}


/**
 * @brief This method performs one transport sweep of all azimuthal angles, 
 *        tracks, segments, polar angles and energy groups.
 * @details The method integrates the flux along each track and updates the 
 *          boundary fluxes for the corresponding output track, while updating 
 *          the scalar flux in each flat source region
 */
void ThreadPrivateSolver::transportSweep() {

    int tid;
    int fsr_id;
    Track* curr_track;
    int azim_index;
    int num_segments;
    segment* curr_segment;    
    segment* segments;
    FP_PRECISION* track_flux;

    log_printf(DEBUG, "Transport sweep with %d OpenMP threads", _num_threads);

    /* Initialize flux in each region to zero */
    flattenFSRFluxes(0.0);

    if (_cmfd->getMesh()->getCmfdOn()) 
      zeroSurfaceCurrents();

    /* Loop over azimuthal angle halfspaces */
    for (int i=0; i < 2; i++) {

        /* Compute the minimum and maximum track IDs corresponding to 
         * this azimuthal angular halfspace */
        int min = i * (_tot_num_tracks / 2);
	int max = (i + 1) * (_tot_num_tracks / 2);
	
	/* Loop over each thread within this azimuthal angle halfspace */
	#pragma omp parallel for private(tid, fsr_id, curr_track, \
	  azim_index, num_segments, segments, curr_segment, \
	  track_flux) schedule(guided)
	for (int track_id=min; track_id < max; track_id++) {

	    tid = omp_get_thread_num();

	    /* Initialize local pointers to important data structures */	
	    curr_track = _tracks[track_id];
	    azim_index = curr_track->getAzimAngleIndex();
	    num_segments = curr_track->getNumSegments();
	    segments = curr_track->getSegments();
	    track_flux = &_boundary_flux(track_id,0,0,0);

	    /* Loop over each segment in forward direction */
	    for (int s=0; s < num_segments; s++) {
	        curr_segment = &segments[s];
		fsr_id = curr_segment->_region_id;
		scalarFluxTally(curr_segment, azim_index, track_flux, 
	                        &_thread_flux(tid,fsr_id,0),true);
	    }

	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, azim_index, true, track_flux);
	    
	    /* Loop over each segment in reverse direction */
	    track_flux += _polar_times_groups;
	    
	    for (int s=num_segments-1; s > -1; s--) {
	        curr_segment = &segments[s];
		fsr_id = curr_segment->_region_id;
		scalarFluxTally(curr_segment, azim_index, track_flux, 
	                        &_thread_flux(tid,fsr_id,0),false);
	    }
	    
	    /* Transfer flux to outgoing track */
	    transferBoundaryFlux(track_id, azim_index, false, track_flux);
	}
    }

    reduceThreadScalarFluxes();
    
    if (_cmfd->getMesh()->getCmfdOn())
      reduceThreadSurfaceCurrents();
    
    return;
}


/**
 * @brief Computes the contribution to the flat source region scalar flux
 *        from a single track segment.
 * @details This method integrates the angular flux for a track segment across
 *        energy groups and polar angles, and tallies it into the flat
 *        source region scalar flux, and updates the track's angular flux.
 * @param curr_segment a pointer to the segment of interest
 * @param track_flux a pointer to the track's angular flux
 * @param fsr_flux a pointer to the temporary flat source region flux buffer
 */
void ThreadPrivateSolver::scalarFluxTally(segment* curr_segment,
					  int azim_index,
					  FP_PRECISION* track_flux,
					  FP_PRECISION* fsr_flux,
					  bool fwd){

    int tid = omp_get_thread_num();
    int fsr_id = curr_segment->_region_id;
    FP_PRECISION length = curr_segment->_length;
    double* sigma_t = curr_segment->_material->getSigmaT();

    /* The average flux along this segment in the flat source region */
    FP_PRECISION deltapsi;
    FP_PRECISION exponential;

    /* Loop over energy groups */
    for (int e=0; e < _num_groups; e++) {

	/* Loop over polar angles */
        for (int p=0; p < _num_polar; p++){
            exponential = computeExponential(sigma_t[e], length, p);
            deltapsi = (track_flux(p,e) - _reduced_source(fsr_id,e)) * exponential;
	    fsr_flux[e] += deltapsi * _polar_weights(azim_index, p);
	    track_flux(p,e) -= deltapsi;
	}
    }

    if (_cmfd->getMesh()->getCmfdOn()){
    	if (curr_segment->_mesh_surface_fwd != -1 && fwd){

    		/* set polar angle * energy group to 0 */
    		int pe = 0;

    		/* loop over energy groups */
    		for (int e = 0; e < _num_groups; e++) {

    			/* loop over polar angles */
    			for (int p = 0; p < _num_polar; p++){

    				/* increment current (polar and azimuthal weighted flux, group)*/
			        _thread_currents(tid,curr_segment->_mesh_surface_fwd,e) += track_flux(p,e)*_polar_weights(azim_index, p)/2.0;

    				pe++;
    			}
    		}
    	}
    	else if (curr_segment->_mesh_surface_bwd != -1 && !fwd){

    		/* set polar angle * energy group to 0 */
    		int pe = 0;

    		/* loop over energy groups */
    		for (int e = 0; e < _num_groups; e++) {

    			/* loop over polar angles */
    			for (int p = 0; p < _num_polar; p++){

    				/* increment current (polar and azimuthal weighted flux, group)*/
			        _thread_currents(tid,curr_segment->_mesh_surface_bwd,e) += track_flux(p,e)*_polar_weights(azim_index, p)/2.0;

    				pe++;
    			}
    		}
    	}
    }

    return;
}


/**
 * @brief Reduces the flat source region scalar fluxes from private thread 
 *        array to a global array.
 */
void ThreadPrivateSolver::reduceThreadScalarFluxes() {

    for (int tid=0; tid < _num_threads; tid++) {
        for (int r=0; r < _num_FSRs; r++) {
            for (int e=0; e < _num_groups; e++) {
                _scalar_flux(r,e) += _thread_flux(tid,r,e);
            }
        }
    }

    return;
}


/**
 * @brief Reduces the surface currents from private thread
 *        array to a global array.
 */
void ThreadPrivateSolver::reduceThreadSurfaceCurrents() {

	for (int tid=0; tid < _num_threads; tid++){
		for (int r=0; r < _num_mesh_cells; r++) {
			for (int s=0; s < 8; s++) {
				for (int e=0; e < _num_groups; e++)
					_surface_currents(r*8+s,e) += _thread_currents(tid,r*8+s,e);
			}
		}
	}

    return;
}
