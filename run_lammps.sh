directories="base_example_T_100 base_example_T_200 base_example_T_300 "
cd ./Examples
for i in $directories; 
do 
cd $i
    lmps -i in.example -l out.example
cd ..
done
cd ..