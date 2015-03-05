package ui;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.HeadlessException;
import java.awt.Image;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.Transparency;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.awt.image.ColorModel;
import java.awt.image.PixelGrabber;
import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JOptionPane;
import javax.swing.JPanel;

import org.apache.commons.collections15.Transformer;
import org.apache.commons.collections15.functors.ConstantTransformer;
import org.apache.commons.collections15.functors.MapTransformer;
import org.apache.commons.collections15.map.LazyMap;
import org.python.modules.synchronize;

import algorithms.communityMining.data.Grouping;
import algorithms.communityMining.topleaders.data.Partitioning;
import edu.uci.ics.jung.algorithms.layout.AggregateLayout;
import edu.uci.ics.jung.algorithms.layout.CircleLayout;
import edu.uci.ics.jung.algorithms.layout.FRLayout;
import edu.uci.ics.jung.algorithms.layout.Layout;
import edu.uci.ics.jung.algorithms.scoring.DegreeScorer;
import edu.uci.ics.jung.algorithms.scoring.util.VertexScoreTransformer;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.SparseMultigraph;
import edu.uci.ics.jung.visualization.GraphZoomScrollPane;
import edu.uci.ics.jung.visualization.VisualizationViewer;
import edu.uci.ics.jung.visualization.control.DefaultModalGraphMouse;
import edu.uci.ics.jung.visualization.decorators.AbstractVertexShapeTransformer;
import edu.uci.ics.jung.visualization.renderers.Renderer;

/**
 * @author rabbanyk
 *
 * @param <V> node type
 * @param <E> edge type
 * 
 * A class for visualizing a network's partitioning/clustering/communities
 *  
 */
public class PartitionView<V,E>  extends JPanel{
	private static final long serialVersionUID = 2475068282663470254L;
	Graph<V, E> graph;
	final Map<E, Double> edgeWeights;
	Grouping<V> groundTruth;
	Grouping<V> partitioning;
	
	
	boolean groupVertices = false;
	int w = 200, h = 200;
	AggregateLayout<V, E> layout;
	VisualizationViewer<V, E> vv;

	@SuppressWarnings("unchecked")
	Map<V, Paint> vertexPaints = LazyMap.<V, Paint> decorate(
			new HashMap<V, Paint>(), new ConstantTransformer(Color.white));
	@SuppressWarnings("unchecked")
	Map<E, Paint> edgePaints = LazyMap.<E, Paint> decorate(
			new HashMap<E, Paint>(), new ConstantTransformer(Color.blue));
	public final Color[] similarColors = { 
			new Color(204,51,51,200),
			new Color(51,204,0,200),
			new Color(0,102,204,200),
			new Color(102,51,102,200),
			new Color(204,204,51,200),
			new Color(51,204,204,200),
			new Color(204,0,204,200),
			new Color(0, 130, 0),
			new Color(255, 193, 193), new Color(28, 141, 255),
			new Color(255, 255, 0), new Color(204,0, 204),
			new Color(255, 130, 5), new Color(102, 51, 0),
			new Color(102, 0, 204), new Color(102, 102, 51),
			new Color(0,0, 204), new Color(153, 0, 0),
			new Color(102, 102, 51), new Color(192, 62, 58) };
//	public final Color[] similarColors = { new Color(0, 130, 0),
//			new Color(255, 193, 193), new Color(28, 141, 255),
//			new Color(255, 255, 0), new Color(204,0, 204),
//			new Color(255, 130, 5), new Color(102, 51, 0),
//			new Color(102, 0, 204), new Color(102, 102, 51),
//			new Color(0,0, 204), new Color(153, 0, 0),
//			new Color(102, 102, 51), new Color(192, 62, 58) };
//	
	
	public PartitionView(Graph<V, E> graph,  Map<E, Double> weights, final Map<V,String> labels,Grouping<V> groundTruth2,
			Grouping<V> partitioning) {
		this.graph = graph;
		this.groundTruth = groundTruth2;
		this.edgeWeights = weights;
		this.partitioning = partitioning;
		this.setBackground(Color.white);
		layout = new AggregateLayout<V, E>(new FRLayout<V, E>(graph));

		vv = new VisualizationViewer<V, E>(layout);
		vv.setBackground(Color.white);
		vv.getRenderer().getVertexLabelRenderer().setPosition(
				Renderer.VertexLabel.Position.CNTR);
		vv.getRenderContext().setVertexLabelTransformer(
				new Transformer<V, String>() {
					public String transform(V arg0) {
						return arg0.toString();// (new Integer(arg0.intValue() +
												// 1)).toString();
					}
				});

		vv.setForeground(Color.black);
		if(labels!=null)
		vv.getRenderContext().setVertexLabelTransformer(
				new Transformer<V, String>() {
					public String transform(V arg0) {
						return labels.get(arg0);// (new Integer(arg0.intValue() +
												// 1)).toString();
					}
				});
		// Tell the renderer to use our own customized color rendering
		vv.getRenderContext().setVertexFillPaintTransformer(
				MapTransformer.<V, Paint> getInstance(vertexPaints));

		vv.getRenderContext().setEdgeDrawPaintTransformer(
				MapTransformer.<E, Paint> getInstance(edgePaints));
		vv.getRenderContext().setEdgeDrawPaintTransformer(new Transformer<E, Paint>() {
   			public Paint transform(E arg0) {
   				return new Color(51,0,0,100);//(arg0.getNormalizedWeight(cst,cet)==0)?null:new Color(51,0,0,100);//Color.black;//new Color(0,0,0);
   			}   	
   			});
		vv.getRenderContext().setVertexStrokeTransformer(new Transformer<V, Stroke>() {
			public Stroke transform(V arg0) {
				return new BasicStroke(0);
			}
		});
		
		if (edgeWeights!=null)
		vv.getRenderContext().setEdgeStrokeTransformer(new Transformer<E, Stroke>() {
			@Override
			public Stroke transform(E input) {
				return new BasicStroke( edgeWeights.get(input).floatValue())  ;
			}
		});
		vv.getRenderContext().setVertexShapeTransformer(
				new VertexShapeSizeAspect(graph,
						new VertexScoreTransformer<V, Integer>(
								new DegreeScorer<V>(graph))));

		DefaultModalGraphMouse<Number, Number> gm = new DefaultModalGraphMouse<Number, Number>();
		vv.setGraphMouse(gm);
        //picking mode with 'p' transferring with 't'
	    vv.addKeyListener(gm.getModeKeyListener());
		
	//	GraphZoomScrollPane graphZoomScrollPane = new GraphZoomScrollPane(vv);
		reDraw();
		this.add(vv);
	}
	private void reDraw() {
		layout.removeAll();
		int i = 0;

		for (V v : vertexPaints.keySet()) {
			vertexPaints.put(v, Color.GRAY);
		}

		if (partitioning!=null)
//			for (V v : partitioning.getOutliers()) {
//				vertexPaints.put(v, Color.WHITE);
//			}

		if (partitioning!=null)
		for (Iterator<Set<V>> cIt = partitioning.getGroups().iterator(); cIt.hasNext();) {

			Set<V> vertices = cIt.next();
			Color c = similarColors[i % similarColors.length];

//			System.err.println(vertices);
			colorCluster(vertices, c);
			if (groupVertices == true) {
				groupCluster(vertices);
			}
			i++;
		}
		for (E e : graph.getEdges()) {
			edgePaints.put(e, Color.black);
		}

		repaint();
	}
		
	private void colorCluster(Set<V> vertices, Color c) {
		for (V v : vertices) {
			vertexPaints.put(v, c);
		}
	}

	private void groupCluster(Set<V> vertices) {
		// System.err.println("____________________________ groupCluster _________________________");

		if (vertices.size() < layout.getGraph().getVertexCount()) {
			Point2D center = layout.transform(vertices.iterator().next());
			Graph<V, E> subGraph = SparseMultigraph.<V, E> getFactory().create();
			for (V v : vertices) {
				subGraph.addVertex(v);
			}
			Layout<V, E> subLayout = new CircleLayout<V, E>(subGraph);
			subLayout.setInitializer(vv.getGraphLayout());
			subLayout.setSize(new Dimension(40, 40));

			layout.put(subLayout, center);
			repaint();
		}
	}
	
	private final class VertexShapeSizeAspect extends	AbstractVertexShapeTransformer<V> implements Transformer<V, Shape> {
		protected boolean scale = true;
		protected Transformer<V, Integer> voltages;
//		protected Graph<V, E> graph;

		public VertexShapeSizeAspect(Graph<V, E> graphIn, Transformer<V, Integer> voltagesIn) {
//			this.graph = graphIn;
			this.voltages = voltagesIn;
			setSizeTransformer(new Transformer<V, Integer>() {

				public Integer transform(V v) {
					if (scale)
						return ((int)((voltages.transform(v) )*1.5 + 10));
					else
						return 1;

				}
			});
		}

		public Shape transform(V arg0) {

			if (groundTruth != null) {
				int s = groundTruth.getGroups().indexOf(groundTruth.getGroup(arg0));
				if (s == -1)
					return factory.getEllipse(arg0);
				if (s < 5) {
					int sides = s + 3;
					return factory.getRegularPolygon(arg0, sides);
				} else
					return factory.getRegularStar(arg0, s);
			} else
				return factory.getEllipse(arg0);

		}

	}

	public void setDividersLocations(int w, int h){
//		mainSplitPane.setDividerLocation((int)(w * .15));
//		graphsSplitPane.setDividerLocation((int)(w * .4));
	}
	
	public synchronized void saveAsFigure(String fileName){
		
		try {
			Image img = createImage(getWidth(), getHeight());
			System.err.println(this.isDisplayable() );
			File saveFile = new File(fileName + ".png");
			Graphics g = img.getGraphics();
			paint(g);
			ImageIO.write(toBufferedImage(img), "png", saveFile);
		//	JOptionPane.showMessageDialog(null,"Image saved to " + saveFile.toString());
			g.dispose();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}
	
	// This method returns a buffered image with the contents of an image
	private BufferedImage toBufferedImage(Image image) {
	    if (image instanceof BufferedImage) {
	        return (BufferedImage)image;
	    }

	    // This code ensures that all the pixels in the image are loaded
	    image = new ImageIcon(image).getImage();

	    // Determine if the image has transparent pixels; for this method's
	    // implementation, see e661 Determining If an Image Has Transparent Pixels
	    boolean hasAlpha = hasAlpha(image);

	    // Create a buffered image with a format that's compatible with the screen
	    BufferedImage bimage = null;
	    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
	    try {
	        // Determine the type of transparency of the new buffered image
	        int transparency = Transparency.OPAQUE;
	        if (hasAlpha) {
	            transparency = Transparency.BITMASK;
	        }

	        // Create the buffered image
	        GraphicsDevice gs = ge.getDefaultScreenDevice();
	        GraphicsConfiguration gc = gs.getDefaultConfiguration();
	        bimage = gc.createCompatibleImage(
	                image.getWidth(null), image.getHeight(null), transparency);
	    } catch (HeadlessException e) {
	        // The system does not have a screen
	    }

	    if (bimage == null) {
	        // Create a buffered image using the default color model
	        int type = BufferedImage.TYPE_INT_RGB;
	        if (hasAlpha) {
	            type = BufferedImage.TYPE_INT_ARGB;
	        }
	        bimage = new BufferedImage(image.getWidth(null), image.getHeight(null), type);
	    }

	    // Copy image to buffered image
	    Graphics g = bimage.createGraphics();

	    // Paint the image onto the buffered image
	    g.drawImage(image, 0, 0, image.getWidth(null), image.getHeight(null), null);
	    g.dispose();

	    return bimage;
	}    
	// This method returns true if the specified image has transparent pixels
	public static boolean hasAlpha(Image image) {
	    // If buffered image, the color model is readily available
	    if (image instanceof BufferedImage) {
	        BufferedImage bimage = (BufferedImage)image;
	        return bimage.getColorModel().hasAlpha();
	    }

	    // Use a pixel grabber to retrieve the image's color model;
	    // grabbing a single pixel is usually sufficient
	     PixelGrabber pg = new PixelGrabber(image, 0, 0, 1, 1, false);
	    try {
	        pg.grabPixels();
	    } catch (InterruptedException e) {
	    }

	    // Get the image's color model
	    ColorModel cm = pg.getColorModel();
	    return cm.hasAlpha();
	}
	
}
