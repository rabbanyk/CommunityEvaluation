
class all_levels {

	
	// this class is used to edit the frames

public:

	all_levels(){};
	~all_levels(){};
	
	
	// set_frames
	void from_partitions_to_frames_no_overlap(int_matrix & short_pt, int_matrix & orig_pt_node);
	int set_frames(deque<int_matrix> & short_parts, deque<int_matrix> & orig_pts_node_no_overlap );
	void original_nodes_in_these_frames(DI & fri, DI & original_n, deque<frame> & last_frames);
	
	// positioning
	
	double minimize_distance_frame(deque<frame*> & Fr, DI & ranking, netvi & frame_network);	
	void frame_location(deque<frame*> & Fr, DI & ranking);
	int swap_homeless_frames(DI & homeless_frames, deque<frame*> & Fr);
	int set_single_frame_positions_no_homeless(deque<frame*> & Fr, int level);
	int set_single_frame_positions(deque<frame*> & Fr, DI & ranking);
	
	int highest_level_set_positions();
	int set_positions(UI );
	void set_radii(frame* A, int level);
	int set_radii_and_positions();
	void cascade(int );
	int reset_frames();
	
	
	// print
	void print_all_frames_inclusions();
	int print_all_frames_radii();
	int print_level_gdf(int level);
	int print_levels_gdf();
	void print_levels_multiple_files();
	void reprint_with_single_nodes(DI & nodes, string filename, map<int, deque<pair<int, double> > > & link_weight_last_level );

	
	// data structures
	netvi original_network;
	deque<deque<frame> > all_level_frames;							// each deque contains all the frames of a certain level
	int_matrix node_mems_across_levels;
	
	
};


#include "all_levels_set_frames.cpp"
#include "all_levels_print.cpp"
#include "all_levels_positioning.cpp"

