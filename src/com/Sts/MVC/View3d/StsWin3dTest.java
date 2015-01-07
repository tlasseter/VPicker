package com.Sts.MVC.View3d;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jul 3, 2008
 * Time: 11:24:40 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsWin3dTest extends JFrame
{
    GLJPanelTest graphicsPanel;
    JPanel objectTreePanelTest;
    JPanel workflowPanelTest;
    JMenuBar menuBar;
    StsCursor3dTest cursor3dTest;
    JPanel statusAreaTest;
    JPanel toolbarPanelTest;

    JSplitPane viewLogSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
    JSplitPane viewTreeSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
    JTabbedPane tabbedPanels = new JTabbedPane();
    JTabbedPane tabbedTextPanels = new JTabbedPane();
    StsJPanelTest logPanel;

    protected Dimension glPanelSize = new Dimension(650, 500); // 3d window
    protected Dimension tabbedPanelSize = new Dimension(270, 500); // object&workflow panels: same height as glPanel
    protected Dimension cursor3dPanelSize = new Dimension(920, 90);
    protected Dimension tabbedTextPanelSize = new Dimension(920, 50); // message panels: width is sum of glPanel and tabbedPanel
    protected Dimension statusAreaSize = new Dimension(920, 25);
    protected Dimension size = new Dimension(920, 665);
    protected Point location = new Point(100, 100);

    public StsNewCursor3dPanel cursor3dPanel = null;

    public StsWin3dTest()
    {
    }

    protected void initialize()
    {
        try
        {
//			setTitle("Version: test");

//			addWindowListener(this); // respond directly to window ops
//            addMouseListener(this);  // used to catch window resize drag events

            graphicsPanel = new GLJPanelTest();
            graphicsPanel.setSize(glPanelSize);
            graphicsPanel.setMinimumSize(new Dimension(0, 0));

            menuBar = new JMenuBar();
//            setJMenuBar(menuBar);

			Container pane = getContentPane();
			pane.setLayout(new BorderLayout());

            statusAreaTest = new JPanel();

            objectTreePanelTest = new JPanel();
            objectTreePanelTest.setSize(tabbedPanelSize);
            workflowPanelTest = new JPanel();
            workflowPanelTest.setSize(tabbedPanelSize);

            tabbedPanels.setPreferredSize(tabbedPanelSize);
            tabbedPanels.setMinimumSize(new Dimension(200, 0));
            tabbedTextPanels.setPreferredSize(tabbedTextPanelSize);
            tabbedTextPanels.setTabPlacement(JTabbedPane.BOTTOM);

            viewTreeSplitPane.setOneTouchExpandable(true);
            viewTreeSplitPane.setContinuousLayout(false);
            viewTreeSplitPane.add(tabbedPanels, JSplitPane.LEFT);
            viewTreeSplitPane.setResizeWeight(0.0); // right component gets all the extra space
            viewTreeSplitPane.setDividerLocation(0.30);

            viewTreeSplitPane.add(graphicsPanel, JSplitPane.RIGHT);

            cursor3dTest = new StsCursor3dTest();
            cursor3dPanel = new StsNewCursor3dPanel(cursor3dPanelSize, cursor3dTest);

            logPanel = new StsJPanelTest();
            logPanel.gbc.fill = GridBagConstraints.BOTH;
            logPanel.add(cursor3dPanel);
            logPanel.add(tabbedTextPanels);
            logPanel.add(statusAreaTest);

            viewLogSplitPane.addPropertyChangeListener(new SplitChangeListener());
            viewLogSplitPane.setOneTouchExpandable(true);
            viewLogSplitPane.setContinuousLayout(false);
            viewLogSplitPane.add(viewTreeSplitPane, JSplitPane.TOP);
            viewLogSplitPane.add(logPanel, JSplitPane.BOTTOM);
            viewLogSplitPane.setResizeWeight(1.0); // top component gets all the extra space

            tabbedPanels.add("Workflow", workflowPanelTest);
            tabbedPanels.add("Objects", objectTreePanelTest);

            toolbarPanelTest = new JPanel();
       //     pane.add(toolbarPanelTest, BorderLayout.NORTH);
       //     pane.add(viewLogSplitPane, BorderLayout.CENTER);
            add(toolbarPanelTest, BorderLayout.NORTH);
            add(viewLogSplitPane, BorderLayout.CENTER);
			pack();
//		    setSize(size);
//		    centerComponentOnScreen(location);
            setVisible(true);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    /*
          protected void addWindowListener(StsWin3dTest window_)
           {
               final StsWin3dTest window = window_;
               WindowListener windowListener = new WindowAdapter()
               {
                   public void windowClosing(WindowEvent e)
                   {
                   }
               };
               window.addWindowListener(windowListener);
           }
    */
    protected void addMouseListener(StsWin3dTest window_)
    {
        final StsWin3dTest window = window_;
        MouseListener mouseListener = new MouseAdapter()
        {
            public void mousePressed(MouseEvent e)
            {
                System.out.println("window: " + getName() + " mouse pressed.");
            }

            public void mouseReleased(MouseEvent e)
            {
                System.out.println("window: " + getName() + " mouse released.");
            }

            public void mouseDragged(MouseEvent e)
            {
                System.out.println("window: " + getName() + " mouse dragged.");
            }
        };
        window.getRootPane().addMouseListener(mouseListener);
    }

    // keep the menu item in sync with the tabs
    class SplitChangeListener implements PropertyChangeListener
    {
        public void propertyChange(PropertyChangeEvent e)
        {
            if (e.getPropertyName().equals(JSplitPane.LAST_DIVIDER_LOCATION_PROPERTY))
            {
                setDividerState();
            }
        }
    }

    private void setDividerState()
    {
        int maxLoc = viewTreeSplitPane.getMaximumDividerLocation();
        boolean selected = viewTreeSplitPane.getDividerLocation() <= maxLoc;
    }

    public void startWindow()
    {

    }


    public StsCursor3dTest getCursor3dTest()
    {
        return cursor3dTest;
    }


    public void display3dCursorDialog(boolean selected)
    {
    }

    public void validate()
    {
    }

    class GLJPanelTest extends GLJPanel implements GLEventListener
    {
        GLJPanelTest()
        {

        }

        public void init(GLAutoDrawable drawable)
        {

        }

        public void display(GLAutoDrawable drawable)
        {

        }

        public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height)
        {

        }

        public void displayChanged(GLAutoDrawable drawable, boolean modeChanged, boolean deviceChanged)
        {

        }
    }

    static public void main(String[] args)
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                StsWin3dTest win3d = new StsWin3dTest();
                win3d.initialize();
//                StsToolkit.createDialog(win3d);
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
    }

}
