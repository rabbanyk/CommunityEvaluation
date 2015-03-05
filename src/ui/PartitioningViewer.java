package ui;

import io.graph.GraphInputStream;
import io.graph.gml.GMLGraphReader;
import io.graph.pajek.PajekGraphReader;
import io.group.GroupingReader;
import io.group.ListGrouingReader;
import io.group.PairGrouingReader;

import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Set;
import java.util.Vector;

import javax.swing.JApplet;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JTabbedPane;
import javax.swing.KeyStroke;

import measure.graph.criteria.Modularity;

import org.apache.commons.collections15.Transformer;

import algorithms.communityMining.data.Grouping;
import algorithms.communityMining.topleaders.TopLeaders;
import algorithms.communityMining.topleaders.data.Partitioning;
import edu.uci.ics.jung.graph.Graph;

/**
 * @author rabbanyk
 * 
 *	A simple Visualization Frame for iTopLeaders Performance
 *	TODO: add other methods as well 
 */
@SuppressWarnings("serial")
public class PartitioningViewer extends JApplet{

	Graph<Integer, Integer> graph;
	Grouping<Integer> groundTruth;
	JTabbedPane tabbedPane;
	
	public PartitioningViewer() {
		setUpMenu();
		setUpBody();
	}
	
	private static StatusBar statusBar ;
	
	public void setUpMenu(){
		
		JMenuBar menuBar;
		JMenu menu;
		JMenuItem menuItem;

		//Create the menu bar.
		menuBar = new JMenuBar();
		this.setJMenuBar(menuBar);

		//Build the first menu.
		menu = new JMenu("File");
		menu.setMnemonic(KeyEvent.VK_F);
		menu.getAccessibleContext().setAccessibleDescription("Menues related to import/export data");
		menuBar.add(menu);
		
		//File Menu
		menuItem = new JMenuItem("Open", KeyEvent.VK_O);
		menuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, ActionEvent.CTRL_MASK));
		menuItem.getAccessibleContext().setAccessibleDescription("importing data");
		
		menuItem.addActionListener(new ActionListener(){
			//Create a file chooser
			final JFileChooser fc = new JFileChooser(new File("./"));

			public void actionPerformed(ActionEvent e) {
			        int returnVal = fc.showOpenDialog(null);

			        if (returnVal == JFileChooser.APPROVE_OPTION) {
			            File file = fc.getSelectedFile();
			            try {
							loadDataSet(file);
						} catch (FileNotFoundException e1) {
							e1.printStackTrace();
						} catch (IOException e1) {
							e1.printStackTrace();
						}
			        } 
			   }

			});
		menu.add(menuItem);
		
		//Build second menu in the menu bar.
		menu = new JMenu("Partition");
		menu.setMnemonic(KeyEvent.VK_P);
		menu.getAccessibleContext().setAccessibleDescription("Partitioning the graph");
		menuBar.add(menu);
	
		menuItem = new JMenuItem("local TopLeader", KeyEvent.VK_L);
		menuItem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_L, ActionEvent.CTRL_MASK));
		menuItem.getAccessibleContext().setAccessibleDescription("performing  top leader data");
		
		menuItem.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
			       showLocalTopLeader();
			   }

			});
		menu.add(menuItem);
		
		menu.addSeparator();	
		
		//this.setJMenuBar(menuBar);

	}
	
	public void setUpBody(){
		tabbedPane = new JTabbedPane();
		this.getContentPane().add(tabbedPane);
		
		statusBar = new StatusBar();
		getContentPane().add(statusBar, java.awt.BorderLayout.SOUTH);
		setStatus("Load the graph by pressing Ctr+o");
	}
	
	//TODO : change to show gt & add this 
	void showLocalTopLeader(){
		System.err.println("Showing partition View");
		
//		LocalTopLeader<Integer, Integer> localLeader = new LocalTopLeader<Integer, Integer>();
		Vector<Vector<Set<Integer>>> vector = new Vector<Vector<Set<Integer>>>();
//		TopLeaders<Integer, Integer> topLeaders = new TopLeaders<Integer, Integer>();
		int max = 0;
		Modularity<Integer, Integer> modularity = new Modularity<Integer, Integer> (graph);
		for (int i = 0; i < 4; i++) {
			vector.add(( new TopLeaders<Integer, Integer>(i+2)).findCommunities(graph, null).getGroups());
			if(modularity.evaluate(vector.get(i)) >= modularity.evaluate(vector.get(max))) max = i;
		}
		//Partitioning< Integer> partitioning = groundTruth;//localLeader.transform(graph);
		Grouping<Integer> partitioning = new Grouping<Integer>();//vector.get(max));//localLeader.transform(graph);
		
		System.err.println(partitioning);
		PartitionView<Integer, Integer> partitionView = new PartitionView<Integer, Integer>(graph, null,null,groundTruth, partitioning);

        tabbedPane.addTab("Partitioned network", partitionView);
        partitionView.setDividersLocations(tabbedPane.getWidth(), tabbedPane.getHeight());
        tabbedPane.setSelectedComponent(partitionView);
	}
	
	public void show(PartitionView<Integer, Integer> partitionView ){
		  tabbedPane.addTab("Partitioned network", partitionView);
	        partitionView.setDividersLocations(tabbedPane.getWidth(), tabbedPane.getHeight());
	        tabbedPane.setSelectedComponent(partitionView);
	}
	
	void loadDataSet(File file) throws FileNotFoundException, IOException{
		
		final GraphInputStream<Integer, Integer> graphReader ;
		
		if(file.getName().contains(".pajak")) graphReader =  new PajekGraphReader<Integer, Integer>();
		else if(file.getName().contains(".gml")) graphReader =  new GMLGraphReader<Integer, Integer>();
		//TODO: correct this
		//else if(file.getName().contains(".pairs") || file.getName().contains(".dat")) graphReader =  new PairsGraphReader<Integer, Integer>(); 
		else  graphReader =  null;

		graph = graphReader.readGraph(new FileInputStream(file), null, null);
		
		System.out.print(graph);
		
		Transformer<String,Integer> vertexTransformer = new Transformer<String, Integer>() {
			public Integer transform(String vertex) {
				return graphReader.transform(vertex);
			}
		};
		
		//Load the ground truth exists
		try{
			if(groundTruth==null){
				GroupingReader<Integer> partitioningReader = new ListGrouingReader<Integer>();
				File gtFile = new File(file.getParent(), file.getName().substring(0,file.getName().indexOf('.')) + ".gt");

				if (gtFile.exists()){ 
					partitioningReader = new ListGrouingReader<Integer>();
				}else { 
					gtFile = new File(file.getParent(),file.getName().replaceAll("network", "community"));
					System.err.println(gtFile.getName());
					partitioningReader = new PairGrouingReader<Integer>();
				}  
				
				groundTruth = partitioningReader.readPartitioning(new FileInputStream(gtFile), vertexTransformer);
				
			}
		}catch (FileNotFoundException e) {
			System.err.println("No ground truth finded for the given data set in:" + file.getParent());
			e.printStackTrace();
		}
	}

	public static void setStatus(String s) {
		statusBar.setMessage(s);
		statusBar.repaint();
	}

	public class StatusBar extends JLabel {
		/** Creates a new instance of StatusBar */
		public StatusBar() {
			super();
			setPreferredSize(new Dimension(100, 16));
			setMessage("Ready");
		}
		public void setMessage(String message) {
			setText(" " + message);
		}
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		PartitioningViewer cd = new PartitioningViewer();
		cd.start();
		JFrame jf = new JFrame();
		jf.getContentPane().add(cd);
		jf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		jf.setLocation((int)screenSize.getWidth()/6, (int)screenSize.getHeight()/6);
		jf.setSize((int)screenSize.getWidth()*2/3,(int)screenSize.getHeight()*2/3);
		jf.pack();
		jf.setVisible(true);
	}

}
