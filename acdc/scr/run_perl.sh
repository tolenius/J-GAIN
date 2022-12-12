#!/bin/bash

if [ "x$1" != "x" ]; then 
  set -e
  cd $1 
  set +e
fi 


set -x 

vapors=("A" "N")

wrapper_output=driver_acdc_J_wrapper.f90
make_output=my_Makefile
make_incl=acdc_include.mk

l_incl_ions=1
l_const_vapor=1


temp_values=(280 285 290)  ## K   will not be used if variable_temperature=1
rh_values=(20 30 40)     ## g/kg  will not be used if variable_temperature=1

variable_temperature=1 # 1-is yes 0 -no 


# Create the Perl option string
perl_opt_0=""
cluster_file="input_"
vapor_suffix="_"
for vapor in "${vapors[@]}"; do
    perl_opt_0+=" --cs_only 1$vapor,0"
    [ $l_const_vapor -eq 1 ] && perl_opt_0+=" --no_eq 1$vapor"
    
    cluster_file+="$vapor"
    vapor_suffix+="$vapor"
done
cluster_file+="narrow_neutral"
if [ $l_incl_ions -eq 1 ]; then
    perl_opt_0+=" --variable_ion_source"
    cluster_file+="_neg_pos"
    vapor_suffix+="_ions"
else
    vapor_suffix+="_noions"
fi
cluster_file+=".inp"



if [ $variable_temperature -eq 1 ]; then 
    ###### generate for none RH and variable temperature 
    perl_opt=$perl_opt_0
	perl_opt+=" --variable_temp"
    sub_suffix="_no_rh" 	
    suffix="_no_rh" 
    # Generate the equations
    append_to_file_names=$vapor_suffix$suffix 
    #perl acdc_2021_02_22_daniel.pl --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --e ./Perl_input/HS298.15K_426clusters2016Apr25.txt --dip ./Perl_input/dip_pol_298.15K_426clusters2016Apr25.txt `echo "$perl_opt --i ./Perl_input/$cluster_file --append_to_file_names $append_to_file_names --append_to_subroutine_names $sub_suffix"`
    perl acdc_2021_09_28.pl --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --e ./Perl_input/HS298.15K_426clusters2016Apr25.txt --dip ./Perl_input/dip_pol_298.15K_426clusters2016Apr25.txt `echo "$perl_opt --i ./Perl_input/$cluster_file --append_to_eq $append_to_file_names --append_to_eq_names $sub_suffix  --no_ambient"`

else 

for temperature in ${temp_values[@]}; do 
  for rh in ${rh_values[@]}; do 

    perl_opt=$perl_opt_0
    
    # Done with input
    ##################################################################################################
    
    
    perl_opt+=" --temperature $temperature"
    perl_opt+=" --rh $rh"
    suffix_rhtemp=$((1000*$temperature + $rh))
    suffix="_$suffix_rhtemp"
    sub_suffix="_$suffix_rhtemp"	
    
    # Generate the equations
    append_to_file_names=$vapor_suffix$suffix 
    #perl acdc_2021_02_22_daniel.pl --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --e ./Perl_input/HS298.15K_426clusters2016Apr25.txt --dip ./Perl_input/dip_pol_298.15K_426clusters2016Apr25.txt `echo "$perl_opt --i ./Perl_input/$cluster_file --append_to_file_names $append_to_file_names --append_to_subroutine_names $sub_suffix"`
    perl acdc_2021_09_28.pl --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --e ./Perl_input/HS298.15K_426clusters2016Apr25.txt --dip ./Perl_input/dip_pol_298.15K_426clusters2016Apr25.txt `echo "$perl_opt --i ./Perl_input/$cluster_file  --append_to_eq $append_to_file_names --append_to_eq_names $sub_suffix  --no_ambient"`
 
  done 
done 


fi 


	
	
	
####################  generate wrapper 


cat > $wrapper_output << EOF
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
EOF





if [ $variable_temperature -eq 1 ]; then 
suffix="_no_rh" 
cat >> $wrapper_output << EOF
    
	external feval$suffix,jeval$suffix,formation$suffix 	
    !
	if (rh < 0.) then 
       call acdc_driver(feval$suffix,jeval$suffix,formation$suffix,&
         neqn,nclust,nout_all,c,cs_ref,&
         temperature,ipr,t_max,t_tot,t_iter,ipar,ok,j_out)
		return 
	else 
	   write(*,*) "positive rh cannot be passed for variable temperature configuration"
	   stop 1
	end if 
EOF

else 

for rh in ${rh_values[@]}; do 
  for temp in ${temp_values[@]}; do 
  suffix=$((1000*$temp + $rh))
cat >> $wrapper_output << EOF
    external feval_$suffix,jeval_$suffix,formation_$suffix 
EOF
  done   
done 


cat >> $wrapper_output << EOF
    !
	combined = int(int(temperature)*1000+int(rh))
	
EOF


for rh in ${rh_values[@]}; do 
  for temp in ${temp_values[@]}; do 
  suffix=$((1000*$temp + $rh))
cat >> $wrapper_output << EOF
	if(combined .eq. $suffix) then 
       call acdc_driver(feval_$suffix,jeval_$suffix,formation_$suffix,&
         neqn,nclust,nout_all,c,cs_ref,&
         temperature,ipr,t_max,t_tot,t_iter,ipar,ok,j_out)
       return
	end if
EOF
  done   
done 


fi 

  
cat >> $wrapper_output << EOF
    write(*,*) ' Not found Relative humidity and Temperature', rh , temperature 
    stop 1
end subroutine acdc_driver_wrapper

EOF


######## Generate make file 

if [ $variable_temperature -eq 1 ]; then 
  suffix="_no_rh" 
  oarr=(acdc_equations${vapor_suffix}$suffix.o)
  farr=(acdc_equations${vapor_suffix}$suffix.f90)
else
   oarr=()
   farr=()
   for rh in ${rh_values[@]}; do 
     for temp in ${temp_values[@]}; do 
       suffix=$((1000*$temp + $rh))
       oarr+=(acdc_equations${vapor_suffix}_$suffix.o) 
	   farr+=(acdc_equations${vapor_suffix}_$suffix.f90)
     done   
   done 
fi 

echo '
FC = gfortran
FCFLAGS = -O0 -g -ffree-line-length-none -cpp -fcheck=bounds -finit-local-zero

run = run_acdc_J_example.f90
get_J = get_acdc_J.f90
driver = driver_acdc_J.f90
driver_wrapper = driver_acdc_J_wrapper.f90
system = acdc_system.f90

' > $make_output
echo " 
run: run_acdc_J.o get_acdc_J.o driver.o driver_wrapper.o acdc_system.o acdc_simulation_setup.o solution_settings.o vode.o vodea.o ${oarr[@]}" >> $make_output
echo "ACDC_OBJ = get_acdc_J.o driver.o driver_wrapper.o acdc_system.o acdc_simulation_setup.o solution_settings.o vode.o vodea.o ${oarr[@]}" > $make_incl



echo '	$(FC) $(FCFLAGS) $^ -o $@

run_acdc_J.o: $(run) get_acdc_J.o acdc_simulation_setup.o
	$(FC) $(FCFLAGS) -c $< -o $@

get_acdc_J.o: $(get_J) driver.o driver_wrapper.o acdc_system.o acdc_simulation_setup.o
	$(FC) $(FCFLAGS) -c $< -o $@

driver.o: $(driver) acdc_system.o acdc_simulation_setup.o solution_settings.o
	$(FC) $(FCFLAGS) -c $< -o $@
	
driver_wrapper.o: $(driver_wrapper) driver.o acdc_system.o acdc_simulation_setup.o solution_settings.o
	$(FC) $(FCFLAGS) -c $< -o $@
	
acdc_system.o: $(system)
	$(FC) $(FCFLAGS) -c $< -o $@
' >> $make_output

if [ $variable_temperature -eq 1 ]; then 
       suffix="_no_rh"
       echo "acdc_equations${vapor_suffix}$suffix.o: acdc_equations${vapor_suffix}$suffix.f90 acdc_simulation_setup.o" >> $make_output
       echo '	$(FC) $(FCFLAGS) -c $< -o $@' >> $make_output
       echo "" >> $make_output
else 
   for rh in ${rh_values[@]}; do 
     for temp in ${temp_values[@]}; do 
       suffix=$((1000*$temp + $rh))
       echo "acdc_equations${vapor_suffix}_$suffix.o: acdc_equations${vapor_suffix}_$suffix.f90 acdc_simulation_setup.o" >> $make_output
       echo '	$(FC) $(FCFLAGS) -c $< -o $@' >> $make_output
       echo "" >> $make_output
   
     done   
   done 
fi 

echo 'acdc_simulation_setup.o: acdc_simulation_setup.f90 acdc_system.o
	$(FC) $(FCFLAGS) -c $< -o $@

solution_settings.o: solvers/solution_settings.f90
	$(FC) $(FCFLAGS) -c $< -o $@

vode.o: solvers/vode.f
	$(FC) $(FCFLAGS) -std=legacy -c $< -o $@

vodea.o: solvers/vodea.f
	$(FC) $(FCFLAGS) -std=legacy -c $< -o $@

.PHONY: clean

clean:
	rm -f *.o *.mod run

cleanall:
	rm -f *.o *.mod run  $(driver_wrapper) $(system) \' >> $make_output
echo "${farr[@]}" >> $make_output
