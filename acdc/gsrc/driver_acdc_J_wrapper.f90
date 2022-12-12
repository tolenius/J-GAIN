subroutine acdc_driver_wrapper(neqn,nclust,nout_all,c,cs_ref,&
    temperature,rh,ipr,t_max,t_tot,t_iter,ipar,ok,j_out)
	use acdc_system, only : n_charges
    use driver_acdc_J, only : acdc_driver 
	implicit none
	
	! Input and output
	! Cluster distribution
	integer, intent(in) :: neqn, nclust, nout_all(n_charges)	! total length of the c vector, the number of actual clusters and indices for outgoing flux
	real(kind(1.d0)), intent(inout) :: c(neqn)			! initial concentrations -> final concentrations (1/m^3)
	! Ambient conditions
	real(kind(1.d0)), intent(in) :: cs_ref				! reference coagulation sink (1/s)
	real(kind(1.d0)), intent(in) :: temperature			! temperature (K)
	real(kind(1.d0)), intent(in) :: rh			        ! relative humidity (g/kg) ! negative for none rh and varable temperature
	real(kind(1.d0)), intent(in) :: ipr					! ion production rate (1/s/m^3)
	! Simulation settings and outcome
	real(kind(1.d0)), intent(in) :: t_max, t_tot		! simulation time and total accumulated time (s)
	real(kind(1.d0)), intent(inout) :: t_iter       	! iteration time step for the Euler method (s)
	integer, intent(inout) :: ipar(4)					! parameters for re-calling the monomer settings and rate constants
	logical, intent(out) :: ok							! .false. if integration failed
	real(kind(1.d0)), intent(out) :: j_out(n_charges)	! simulated formation rates (neutral, neg and pos) (1/s/m^3)
	integer :: combined
	!
    
	external feval_no_rh,jeval_no_rh,formation_no_rh 	
    !
	if (rh < 0.) then 
       call acdc_driver(feval_no_rh,jeval_no_rh,formation_no_rh,&
         neqn,nclust,nout_all,c,cs_ref,&
         temperature,ipr,t_max,t_tot,t_iter,ipar,ok,j_out)
		return 
	else 
	   write(*,*) "positive rh cannot be passed for variable temperature configuration"
	   stop 1
	end if 
    write(*,*) ' Not found Relative humidity and Temperature', rh , temperature 
    stop 1
end subroutine acdc_driver_wrapper

