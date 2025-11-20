directories="PuO2_T_100 PuO2_T_200 "
cd .
for i in $directories; 
do 
cd $i
    srun --distribution=block:block --hint=nomultithread/work/e89/e89/wd324/Programs/lammps/build/lmp -i in.PuO2 -l out.PuO2
cd ..
done
cd ..