#!/bin/bash

set -e # stop when an error occurs
set -x # print the executed code

vapors=("A" "N")
l_const_vapor=1
l_incl_ions=1
variable_temperature=1 # 1-is yes 0 -no 

temp_values=(280 285 290)  ## (K)   will not be used if variable_temperature=1
rh_values=(20 30 40)       ## (%)  will not be used if variable_temperature=1

## Files below should be in the input Folder
cluster_file_suffix="_input" # Make sure that the cluster file is compatible with the above choices.
dip_file="dip_pol_298.15K_426clusters2016Apr25.txt" # Uded only if l_incl_ions is 1
hs_file="HS298.15K_426clusters2016Apr25.txt"


######################## User should not make changes below ########################

SCRIPT_PATH=$(realpath $0)
ACDC_DIR=$(dirname $SCRIPT_PATH)


BIN_AD=$ACDC_DIR/build/bin
LIB_AD=$ACDC_DIR/build/bin
INCLUDE_AD=$ACDC_DIR/build/include
WORK_AD=$ACDC_DIR/build/work



##
# A: Absolute
# R: Relative
# D: Directory/Folder
# F: File


SCR_AD=$ACDC_DIR/scr
SRC_AD=$ACDC_DIR/src
GSRC_AD=$ACDC_DIR/gsrc
INPUT_AD=$ACDC_DIR/input

ACDC_PERL_AF=$SCR_AD/acdc_2022_11_24.pl
DIP_FILE_AF=$INPUT_AD/$dip_file
HS_FILE_AF=$INPUT_AD/$hs_file

##################################################


export BIN_AD=$ACDC_DIR/build/bin
export LIB_AD=$ACDC_DIR/build/lib
export INCLUDE_AD=$ACDC_DIR/build/include
export WORK_AD=$ACDC_DIR/build/work




[ -d $GSRC_AD ] && rm -rf $GSRC_AD
[ -d $BIN_AD ] && rm -rf $BIN_AD
[ -d $LIB_AD ] && rm -rf $LIB_AD
[ -d $INCLUDE_AD ] && rm -rf $INCLUDE_AD
[ -d $WORK_AD ] && rm -rf $WORK_AD

mkdir -p $GSRC_AD
mkdir -p $BIN_AD
mkdir -p $LIB_AD
mkdir -p $INCLUDE_AD
mkdir -p $WORK_AD

#########################################################################################################
pushd $GSRC_AD
# Create the Perl option string
perl_opt_0=""
vapor_suffix="_"


wrapper_output=driver_acdc_J_wrapper.f90



# Create the Perl option string
perl_opt_0=""
cluster_file=""
for vapor in "${vapors[@]}"; do
    perl_opt_0+=" --cs_only 1$vapor,0"
    [ $l_const_vapor -eq 1 ] && perl_opt_0+=" --no_eq 1$vapor"
    
    cluster_file+="$vapor"
    vapor_suffix+="$vapor"
done
cluster_file+="_neutral"
dip_option=""
if [ $l_incl_ions -eq 1 ]; then
    perl_opt_0+=" --variable_ion_source"
    cluster_file+="_neg_pos"
    vapor_suffix+="_ions"
	dip_option="--dip $DIP_FILE_AF"
else
    vapor_suffix+="_noions"
fi


cluster_file+=$cluster_file_suffix".inp"


 

CLUSTER_FILE_AF=$INPUT_AD/$cluster_file


if [ $variable_temperature -eq 1 ]; then 
    ###### generate for none RH and variable temperature 
    perl_opt=$perl_opt_0
	perl_opt+=" --variable_temp"
    sub_suffix="_no_rh" 	
    suffix="_no_rh" 
    # Generate the equations
    append_to_file_names=$vapor_suffix$suffix 
    perl $ACDC_PERL_AF --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --exp_loss_ref_size 0.55 --e $HS_FILE_AF $dip_option `echo "$perl_opt --i $CLUSTER_FILE_AF --append_to_eq $append_to_file_names --append_to_eq_names $sub_suffix  --no_ambient"`

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
    perl $ACDC_PERL_AF --fortran --save_outgoing --variable_cs --cs exp_loss --exp_loss_exponent -1.6 --exp_loss_ref_size 0.55 --e $HS_FILE_AF $dip_option `echo "$perl_opt --i $CLUSTER_FILE_AF  --append_to_eq $append_to_file_names --append_to_eq_names $sub_suffix  --no_ambient"`
 
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
	real(kind(1.d0)), intent(in) :: rh			        ! relative humidity (%) ! negative for none rh and varable temperature
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

popd 

############## Compile 


if [ $variable_temperature -eq 1 ]; then 
  suffix="_no_rh" 
  farr=(acdc_equations${vapor_suffix}_$suffix.f90)
else
   farr=()
   for rh in ${rh_values[@]}; do 
     for temp in ${temp_values[@]}; do 
       suffix=$((1000*$temp + $rh))
	   farr+=(acdc_equations${vapor_suffix}_$suffix.f90)
     done   
   done 
fi 



FC="gfortran"
FCFLAGS="-O0 -g -static -ffree-line-length-none -cpp -fcheck=bounds -finit-local-zero"


#
F=solution_settings
$FC $FCFLAGS -c $SRC_AD/solvers/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD  

F=vode
$FC $FCFLAGS -c $SRC_AD/solvers/$F.f -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD 
 
F=vodea
$FC $FCFLAGS -c $SRC_AD/solvers/$F.f -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD 
 
	   
	  
#
F=acdc_system
$FC $FCFLAGS -c $GSRC_AD/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD  

#
F=acdc_simulation_setup
$FC $FCFLAGS -c $SRC_AD/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD  



if [ $variable_temperature -eq 1 ]; then 
       suffix="_no_rh"
	   F=acdc_equations${vapor_suffix}_$suffix
       $FC $FCFLAGS -c $GSRC_AD/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD    
 
else 
   for rh in ${rh_values[@]}; do 
     for temp in ${temp_values[@]}; do 
       suffix=$((1000*$temp + $rh))
	   F=acdc_equations${vapor_suffix}_$suffix
       $FC $FCFLAGS -c $GSRC_AD/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD    
     done   
   done 
fi 



#
F=driver_acdc_J
$FC $FCFLAGS -c $SRC_AD/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD  

#
F=driver_acdc_J_wrapper
$FC $FCFLAGS -c $GSRC_AD/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD  

#
F=get_acdc_J
$FC $FCFLAGS -c $SRC_AD/$F.f90 -o $WORK_AD/$F.o  -I$INCLUDE_AD  -J$INCLUDE_AD  
	
	
ar cr $LIB_AD/libacdc.a $WORK_AD/*.o 

#
F=run_acdc_J_example
$FC $FCFLAGS   $SRC_AD/$F.f90 -o $BIN_AD/$F.exe -lacdc -I$INCLUDE_AD  -J$INCLUDE_AD  -L$LIB_AD
	

