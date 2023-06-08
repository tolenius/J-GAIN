#!/bin/bash
script_dir=`pwd`

work_dir=$script_dir/work

opt=$1

echo "This script is to run programs in serial mode."
echo "Usage:"
echo "bash run_me.sh  [option]"
echo
echo "Please Run the options in order:"
echo "Select Option:"
echo "0- Clean all."
echo "1- Compile ACDC."
echo "2- Compile the Generator."
echo "3- Compile the Interpolator."
echo "4- Run the Generator Reference Case."
echo "5- Run the Generator Experiment Case."
echo "6- Run the Interpolator."
echo "7- Plot."
echo "8- Run all above in order."


if [ "x$opt" == "x" ]; then 
    echo "Inter one of the options above:"
    read opt 
fi 

acd_dir=$work_dir/acdc
gen_dir=$work_dir/generator
int_dir=$work_dir/interpolator

acd_origin_dir=$script_dir/../acdc
gen_origin_dir=$script_dir/../generator
int_origin_dir=$script_dir/../interpolator



#cp -rf $acd_origin_dir $work_dir/ 
#cp -rf $gen_origin_dir $work_dir/ 
#cp -rf $int_origin_dir $work_dir/ 

 
gen_bin=$gen_dir/build.serial/bin/jgain_serial.exe


out_dir=$work_dir/output
ref_dir=$out_dir/ref
exp_dir=$out_dir/exp
t2t_dir=$out_dir/t2t
t2t_bin=$t2t_dir/t2t.exe

if [ "$opt" == "8" ]; then 
  bash $0 0
  bash $0 1
  bash $0 2
  bash $0 3
  bash $0 4
  bash $0 5
  bash $0 6
  bash $0 7
  exit 0 
fi 
case $opt in
    0)
        rm -rf $work_dir
        ;;
    1)
        [ -d $work_dir ] ||  mkdir -p $work_dir
        rm -rf $acd_dir
        cp -rf $acd_origin_dir $work_dir/
        pushd $acd_dir
        bash build-acdc.sh
        echo "Done!"
        exit 0 
        ;;
    2)
        rm -rf $gen_dir
        cp -rf $gen_origin_dir $work_dir/ 
        pushd $gen_dir
        bash build-gen.sh serial
        echo "Done!"
        exit 0 
        ;;
    3)
        rm -rf $int_dir
        cp -rf $int_origin_dir $work_dir/ 
        pushd $int_dir
        bash build-interpolator.sh 
        echo "Done!"
        exit 0
        ;;
    4)
        
        [ -f $gen_bin ] || echo "Compile the Generator" || exit 1 
        [ -d $ref_dir ] || mkdir -p $ref_dir
        ln -sf $script_dir/nam/namelist.ref.gen $ref_dir/namelist.gen
        pushd $ref_dir
        $gen_bin        
        ;;
    5)
        [ -f $gen_bin ] || echo "Compile the Generator" || exit 1 
        [ -d $exp_dir ] || mkdir -p $exp_dir
        ln -sf $script_dir/nam/namelist.exp.gen $exp_dir/namelist.gen
        pushd $exp_dir
        $gen_bin  
		exit 0
        ;;    
    6)
        rm -rf $t2t_dir
        mkdir -p  $t2t_dir
        echo "Compiling table to table interpolater."
        gfortran  -g -ffree-line-length-none -march=native -fcheck=bounds -I$int_dir/build/include $int_dir/build/obj/*.o $script_dir/src/table_to_table_interpolator.f90 -o  $t2t_bin  


        ln -sf $script_dir/nam/namelist.t2t $t2t_dir/namelist.t2t
        pushd $t2t_dir
        pwd
        $t2t_bin
		# just cleanup to be used by python 
		for i in `seq 1 10`
		do
		   sed -i "s/  / /g" t2t.interp-logj2.jump1  # the same as in the mae list 
		done 
		sed -i "s/ /,/g" t2t.interp-logj2.jump1
		sed -i "1i ,num,H2SO4,NH3,ref,exp," t2t.interp-logj2.jump1
		exit 0 	
        ;;    
        
    7)
        cd $out_dir
        python3 $script_dir/src/plot2d.py 
		exit 0
        ;;        
        
    *)
        echo "Please inter a valid option."
        ;;
esac

