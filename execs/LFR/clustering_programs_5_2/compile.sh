


rm -r bin
mkdir bin

Main_folder=Sources


source_folder=$Main_folder/OSLOM_files
visual_folder=$Main_folder/Cvis_sources_1_5_netcom


echo "Compiling OSLOM undirected (oslom_undir) ..."
echo "g++ -o oslom_undir $source_folder/main_undirected.cpp -O3 -Wall"
g++ -o oslom_undir $source_folder/main_undirected.cpp -O3 -Wall
mv oslom_undir bin/.

echo ""
echo "Compiling OSLOM directed (oslom_dir) ..."
echo "g++ -o oslom_dir $source_folder/main_directed.cpp -O3 -Wall"
g++ -o oslom_dir $source_folder/main_directed.cpp -O3 -Wall
mv oslom_dir bin/.


echo ""
echo "Compiling cvis_undir"
echo "g++ -o cvis_undir $visual_folder/main_undir.cpp -O3 -Wall"
g++ -o cvis_undir $visual_folder/main_undir.cpp -O3 -Wall
mv cvis_undir bin/.


echo ""
echo "Compiling cvis_dir"
echo "g++ -o cvis_dir $visual_folder/main_dir.cpp -O3 -Wall"
g++ -o cvis_dir $visual_folder/main_dir.cpp -O3 -Wall
mv cvis_dir bin/.



echo ""
echo "Compiling infomap_undirected ..."
cd $Main_folder/infomap_undir/
make clean
make
g++ -o infomap_scr infomap_scr.cpp -O3
cd ../..
mv $Main_folder/infomap_undir/infomap infomap_undir
mv $Main_folder/infomap_undir/infomap_scr infomap_undir_script
mv infomap_undir bin/.
mv infomap_undir_script bin/.


echo ""
echo "Compiling infomap_directed ..."
cd $Main_folder/infomap_dir/
make clean
make
g++ -o infomap_scr infomap_scr.cpp -O3
cd ../..
mv $Main_folder/infomap_dir/infomap infomap_dir
mv $Main_folder/infomap_dir/infomap_scr infomap_dir_script
mv infomap_dir bin/.
mv infomap_dir_script bin/.


echo ""
echo "Compiling hier infomap_undirected ..."
cd $Main_folder/infohiermap_undir/
make clean
make
g++ -o infomap_scr infomap_scr.cpp -O3
cd ../..
mv $Main_folder/infohiermap_undir/infohiermap infohiermap_undir
mv $Main_folder/infohiermap_undir/infomap_scr infohiermap_undir_script
mv infohiermap_undir bin/.
mv infohiermap_undir_script bin/.

echo ""
echo "Compiling hier infomap_directed ..."
cd $Main_folder/infohiermap_dir/
make clean
make
g++ -o infomap_scr infomap_scr.cpp -O3
cd ../..
mv $Main_folder/infohiermap_dir/infomap infohiermap_dir
mv $Main_folder/infohiermap_dir/infomap_scr infohiermap_dir_script
mv infohiermap_dir bin/.
mv infohiermap_dir_script bin/.


echo ""
echo "Compiling louvain  method ..."
cd $Main_folder/sig_louvain/
g++ main_undirected.cpp -Wall -O3 -o a.out
cd ../..
mv $Main_folder/sig_louvain/a.out louvain_method
mv louvain_method bin/.


echo ""
echo "Compiling LPM"
cd Sources/reka_et_al/
g++ -o lpm main_undirected.cpp -O3
cd ../..
mv Sources/reka_et_al/lpm .
mv lpm bin/.



echo "compiling consensus program"
g++ -o pr_consensus Sources/consensus_matrix/consensus_pr/from_partitions_to_network.cpp -O3
mv pr_consensus bin/.



echo "compiling modularity optimization (simulated annealing)"
g++ Sources/modopt/main.cpp -O3 -o bin/modopt



#echo ""
#echo "Compiling general program"
#echo "$Main_folder/main.cpp -o select -O3"
#g++ $Main_folder/main.cpp -o select -O3



