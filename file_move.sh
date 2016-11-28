mkdir 1_combined_energies  2_Len_Jones_only  3_C_only
mv *CLJ* 1_combined_energies/
mv *coul* 2_Len_Jones_only/
mv *LJ* 3_C_only/

cd 1_combined_energies/
python /home/lisa/Desktop/1_Project/0_analysis/00_scripts_GIT/sire_edecomp/file_rename.py

cd ../2_Len_Jones_only/
python /home/lisa/Desktop/1_Project/0_analysis/00_scripts_GIT/sire_edecomp/file_rename.py

cd ../3_C_only/
python /home/lisa/Desktop/1_Project/0_analysis/00_scripts_GIT/sire_edecomp/file_rename.py
