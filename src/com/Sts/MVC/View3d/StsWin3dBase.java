package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
  * <p>Description: </p>
  * <p>Copyright: Copyright (c) 2002</p>
  * <p>Company: S2S Systems LLC</p>
  * @author unascribed
  * @version 1.1
  */

import com.Sts.Actions.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Icons.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsWin3dBase extends StsWindow implements StsSerializable, ComponentListener
 {
     public StsGLPanel3d[] glPanel3ds = new StsGLPanel3d[0];
     //    private StsCursor3d cursor3d = null;
     /** -1 indicates this is a singleViewPanel, otherwise index is index of corresponding multiViewPanel ( 0 to nMultiViewPanels-1 */
     public StsMultiViewPanel currentViewPanel;
     public String windowName = new String("Main Window");
     public boolean isLocked = true;
     public boolean isMouseTracking = false;
     public boolean isCursorTracked = false;
     public StsObject[] selectedObjects = null;
     public StsToolbar[] toolbars = new StsToolbar[0];

     protected Dimension viewPanelSize = new Dimension(650, 550); // 3d window
     protected Dimension tabbedPanelSize = new Dimension(270, 550); // object & workflow panels: same height as glPanel
     protected Dimension cursor3dPanelSize = new Dimension(920, 75);
     protected Dimension tabbedTextPanelSize = new Dimension(920, 50); // message panels: width is sum of glPanel and tabbedPanel
     protected Dimension statusAreaSize = new Dimension(920, 25);
     protected Dimension toolbarPanelSize = new Dimension(920, 50);
     protected Dimension size = new Dimension(920, 750);
     // protected Dimension splashPanelSize = new Dimension(920, 550);
     protected Point location = null;

     protected transient StsGLPanel3d glPanel3d = null; // convenience copy of glPanel for singleViewPanel

     transient Container pane;
     transient GridBagConstraints gbc = new GridBagConstraints();

     //    transient StsMultiSplitPane multiSplitPane;

     transient public StsStatusArea statusArea;

     transient boolean captureNextDraw = false;
     transient boolean captureMovie = false;
     transient int captureType = GRAPHIC;

     transient protected JFrame parent;
     transient StsCheckboxMenuItem[] menubarToolbarItems;

     transient protected StsToolbar currentToolbar;
     //    transient protected boolean passed = true;

     transient protected JPanel toolbarPanel = null;

     // private transient StsWindowActionManager actionManager = null;
     transient public StsWin3dFull parentWindow = null;

     /** use view shifts if true */
     transient private boolean viewShiftOn = true;

     transient public static final int LEFT = 0;
     transient public static final int RIGHT = 1;

     static public final int DESKTOP = 0;
     static public final int WINDOW = 1;
     static public final int GRAPHIC = 2;

     static public final Dimension size00 = new Dimension(0, 0);

     public StsWin3dBase()
     {
         parent = this;
     }

     public StsWin3dBase(StsModel model, StsWin3dBase launchWindow)
     {
         this.model = model;
         //        this.actionManager = actionManager;
         parent = this;
         actionManager = launchWindow.getActionManager();
         this.parentWindow = (model.getWindowFamily(launchWindow)).getParent();
         model.addWindow(parentWindow, this);
         addComponentListener(this);
         this.isLocked = launchWindow.isLocked;
         try
         {
             windowName = model.getWindowFamilyTitle(this);
             constructWindow();
             //             initializeWindowLayout();
             setMouseTracking(launchWindow.isMouseTracking);
             // start();
         }
         catch (Exception e)
         {
             e.printStackTrace();
         }
     }

     public StsActionManager getActionManager()
     {
         return actionManager;
     }

     public void setMouseTracking(boolean isMouseTracking)
     {
         this.isMouseTracking = isMouseTracking;
         if (isMouseTracking)
             actionManager.startAction(StsCursor3dReadout.class);
         else
             actionManager.endAction(StsCursor3dReadout.class);
     }

     public boolean start()
     {
         try
         {
             setFocusable(true);
             pack();
             setLocation();
             initAndStartGL();
             setVisible(true);
             return true;
         }
         catch (Exception e)
         {
             e.printStackTrace();
             return false;
         }
     }

     private void initAndStartGL()
     {
		 GraphicsDevice g = getGraphicsDevice();
         for (StsGLPanel3d glPanel3d : glPanel3ds)
             glPanel3d.initAndStartGL(g);
     }

     public void repaint()
     {
         super.repaint();
         if (Main.isDrawDebug) System.out.println("StsWin3dBase.repaint() " + getName());
     }

     /*
         public void validate()
         {
             glPanel3d.reshaped = false;
             super.validate();
         }
     */
     //Component initialization from DB
     public void initializeTransients(StsModel model, StsWin3dFull parentWindow)
     {
         this.model = model;
         this.parentWindow = parentWindow;
         initializeActionManager();
         initializeCursor3dTransients(model, this);
         initializeViewPanel();
         constructWindow();
         initializeWindowLayout();
         // checkAddDisplayableSections();
         // start();
     }

     public void initializeActionManager()
     {
         actionManager = parentWindow.getActionManager();
     }

     protected void initializeCursor3dTransients(StsModel model, StsWin3dBase win3d)
     {
     }

     protected void constructWindow()
     {
         try
         {
             initializeDisplayIcon();
             setTitle(model.getName() + " " + windowName);
             addWindowCloseOperation();
             addMouseListener(this);  // used to catch window resize drag events
             //            initializeBorderLayout();
             checkCreateToolbars();
             //            centerComponentOnScreen();
         }
         catch (Exception e)
         {
             StsException.outputException("StsWin3dBase.classInitialize() failed.",
                 e, StsException.WARNING);
         }
     }

     public void initializeWindowLayout()
     {
         //        initializeMultiViewSplitPane();

         //       initializeBorderLayout();
         addBorderComponent(toolbarPanel, BorderLayout.NORTH);
         addBorderComponent(currentViewPanel, BorderLayout.CENTER);
         // addBorderComponent(multiSplitPane, BorderLayout.CENTER);

         // currentViewPanel.setMinimumSize(new Dimension(0,0));
         currentViewPanel.setPreferredSize(viewPanelSize);
         //        currentViewPanel.setMinimumSize(new Dimension(200, 200));
     }

     /*
         protected void initializeMultiViewSplitPane()
         {
             multiSplitPane = new StsMultiSplitPane();
             multiSplitPane.setLayout(new StsMultiSplitPaneLayout());
             multiSplitPane.setDividerSize(3);
             multiSplitPane.addRootLeafComponent("view", currentViewPanel);
         }
     */
     protected void initializeDisplayIcon()
     {
         if (!System.getProperty("os.name").startsWith("Windows")) return;
         ImageIcon imageIcon = StsIcon.createIcon(Main.vendorName + ".gif");
         if (imageIcon == null) return;
         Image image = imageIcon.getImage();
         setIconImage(image);
     }

     protected void initializeBorderLayout()
     {
         //          pane = getContentPane();
         //          pane.setLayout(new BorderLayout());
     }

     protected void addBorderComponent(Component comp, Object position)
     {
         add(comp, position);
         //          pane.add(comp, position);
     }

     protected void initializeGridBagLayout()
     {
         pane = getContentPane();
         pane.setLayout(new GridBagLayout());
         gbc = new GridBagConstraints();
         gbc.gridx = 0;
         gbc.gridy = 0;
         //       gbc.gridwidth = 1;
         //       gbc.gridheight = 1;
         gbc.weightx = 1.0;
         gbc.weighty = 1.0;
         //        gbc.anchor = GridBagConstraints.CENTER;
         gbc.fill = GridBagConstraints.BOTH;
     }

     protected void addGridBagRowComponent(Component comp, double weighty)
     {
         gbc.weighty = weighty;
         pane.add(comp, gbc);
         gbc.gridx = 0;
         gbc.gridy += 1;
         gbc.weighty = 1.0;
     }

     protected void addGridBagRowComponent(Component comp)
     {
         pane.add(comp, gbc);
         gbc.gridx = 0;
         gbc.gridy += 1;
     }

     protected void setLocation()
     {
         if (location == null)
             StsToolkit.centerComponentOnScreen(this); // Center the window
         else
             setLocation(location);
     }

     protected void checkCreateToolbars()
     {
         toolbarPanel = new JPanel(new StsFlowLayout(FlowLayout.LEFT));
         if (toolbars.length == 0)
         {
             createAddToolbars();
             getModel().workflowPlugIn.addAdditionalToolbars(this);
         }
         else
         {
             initializeToolbars();
             addToolbars();
         }
         addToolbarMouseListener();
         createToolbarMenuItems();
     }

     protected void createAddToolbars()
     {
         addToolbar(new StsDefaultToolbar(this));
         addToolbar(new StsViewSelectToolbar(this));
         addToolbar(new StsMouseActionToolbar(this));
         addToolbar(new StsTimeActionToolbar(this));
         addToolbar(StsComboBoxToolbar.constructor(this));
     }

     public void checkAddTimeActionToolbar()
     {
         if (hasToolbarNamed(StsTimeActionToolbar.NAME)) return;
         StsTimeActionToolbar timeActionToolbar = new StsTimeActionToolbar(this);
         addToolbar(timeActionToolbar);
     }

     protected void createToolbarMenuItems()
     {
         int nToolbars = toolbars.length;
         menubarToolbarItems = new StsCheckboxMenuItem[nToolbars];
         for (int n = 0; n < nToolbars; n++)
         {
             StsCheckboxMenuItem menuItem = new StsCheckboxMenuItem();
             menuItem.setSelected(toolbars[n].isVisible());
             menubarToolbarItems[n] = menuItem;
             menuItem.setMenuItemListener(toolbars[n].getName(), this, "toolbarSelected");
         }
     }

     public void toolbarSelected(ItemEvent e)
     {
         StsCheckboxMenuItem menuItem = (StsCheckboxMenuItem) e.getItem();
         String toolbarName = menuItem.getText();
         StsToolbar toolbar = StsToolbar.getToolbarNamed(toolbars, toolbarName);
         if (toolbar == null) return;
         boolean visible = (e.getStateChange() == ItemEvent.SELECTED);
         toolbar.setVisible(visible);
     }

     public void setName(String name)
     {
         windowName = name;
     }

     public String getName()
     {
         return windowName;
     }

     public void outputImage(Integer type)
     {
         captureNextDraw = true;
         captureType = type.intValue();
         win3dDisplay();
     }

     public void outputMovie()
     {
         if (!captureMovie)
             captureMovie = true;
         else
             captureMovie = false;
     }

     public boolean captureMovie() { return captureMovie; }

     public boolean captureNextDraw() { return captureNextDraw; }

     public void setCaptureNextDraw(boolean value) { captureNextDraw = value; }

     public int getCaptureType() { return captureType; }

     protected void initializeToolbars()
     {
  //       hackCheckAddComboBoxToolbar();
         for (int n = 0; n < toolbars.length; n++)
         {
             toolbars[n].initialize(this);
         }
     }

     private void hackCheckAddComboBoxToolbar()
     {
         StsToolbar toolbar = StsToolbar.getToolbarNamed(toolbars, StsComboBoxToolbar.NAME);
         if(toolbar != null) return;
         toolbar = StsComboBoxToolbar.constructor(this);
         toolbars = (StsToolbar[]) StsMath.arrayAddElement(toolbars, toolbar);
         addToolbar(toolbar);
     }

     protected void addToolbarMouseListener()
     {
         toolbarPanel.addMouseListener
             (
                 new MouseAdapter()
                 {
                     public void mousePressed(MouseEvent e)
                     {
                         mousePressedInToolbar(e);
                     }
                 }
             );
     }

     protected void mousePressedInToolbar(MouseEvent e)
     {
         Object source = e.getSource();
         if (source != toolbarPanel) return;

         int mods = e.getModifiers();
         if ((((mods & InputEvent.BUTTON3_MASK) != 0) && (e.isShiftDown())) || ((mods & InputEvent.BUTTON2_MASK) != 0))
         {
             JPopupMenu popupMenu = getToolbarPopupMenu();
             toolbarPanel.add(popupMenu);
             popupMenu.show((Component) toolbarPanel, e.getX(), e.getY());
         }
     }

     public JPopupMenu getToolbarPopupMenu()
     {
         JPopupMenu popupMenu = new JPopupMenu("Toolbars");
         for (int n = 0; n < menubarToolbarItems.length; n++)
             popupMenu.add(menubarToolbarItems[n].copy());
         return popupMenu;
     }

     public void checkAddToolbar(StsToolbar toolbar)
     {
         if (!hasToolbarNamed(toolbar.getName()))
             addToolbar(toolbar);
         toolbar.setVisible(true);
         setToolbarMenuState(toolbar, true);
         repaintToolbar();
     }

     public StsToolbar getToolbarNamed(String name)
     {
         return StsToolbar.getToolbarNamed(toolbars, name);
     }
/*
    public StsToolbar getToolbarNamed(String name)
     {
         StsToolbar toolbar = StsToolbar.getToolbarNamed(toolbars, name);
         if(toolbar != null) return toolbar;
         if(name != StsComboBoxToolbar.NAME) return null;
         toolbar = StsComboBoxToolbar.constructor(this);
         toolbars = (StsToolbar[]) StsMath.arrayAddElement(toolbars, toolbar);
         addToolbar(toolbar);
         return toolbar;
     }
 */
     public StsMouseActionToolbar getMouseActionToolbar()
     {
         return (StsMouseActionToolbar) getToolbarNamed(StsMouseActionToolbar.NAME);
     }

     public StsToolbar[] getToolbars()
     {
         Component[] c = toolbarPanel.getComponents();
         int n = 0;
         for (int i = 0; i < c.length; i++)
         {
             if (c[i] instanceof StsToolbar)
             {
                 n++;
             }
         }
         if (n <= 0) return null;
         StsToolbar[] t = new StsToolbar[n];
         n = 0;
         for (int i = 0; i < c.length; i++)
         {
             if (c[i] instanceof StsToolbar)
             {
                 t[n] = (StsToolbar) c[i];
                 n++;
             }
         }

         return t;
     }

     public StsToolbar validateToolbars()
     {
         Component[] c = toolbarPanel.getComponents();
         for (int i = 0; i < c.length; i++)
         {
             if (c[i] instanceof StsToolbar)
             {
                 StsToolbar tb = (StsToolbar) c[i];
                 tb.validateState();

                 if (tb.getName().equals("Cross Plot Toolbar"))
                 {
                     getModel().getNObjects(StsCrossplot.class);
                     if (getModel().getNObjects(StsCrossplot.class) <= 0)
                         closeToolbar(tb);
                 }

             }
         }
         return null;
     }

     public boolean hasToolbarNamed(String name)
     {
         return getToolbarNamed(name) != null;
     }

     public boolean selectToolbarItem(String toolbarName, String itemName, boolean selected)
     {
         StsToolbar toolbar = getToolbarNamed(toolbarName);
         if (toolbar == null) return false;
         return toolbar.selectToolbarItem(itemName, selected);
     }

     public Component getToolbarComponentNamed(String toolbarName, String componentName)
     {
         StsToolbar toolbar = getToolbarNamed(toolbarName);
         if (toolbar == null) return null;
         return toolbar.getComponentNamed(componentName);
     }

     public void hideToolbarPanel()
     {
         toolbarPanel.setVisible(false);
     }

     public void showToolbarPanel()
     {
         toolbarPanel.setVisible(true);
     }

     public void addToolbar(StsToolbar toolbar)
     {
         if (toolbar == null) return;
         addToolbar(toolbar, LEFT);
     }

     public void addToolbars()
     {
         int nToolbars = toolbars.length;

         for (int n = 0; n < nToolbars; n++)
         {
             if(Main.viewerOnly)
             //{
                toolbars[n].forViewOnly();
                //if(toolbars[n].forViewOnly())
                    //toolbarPanel.add(toolbars[n]);
             //}
             //else
                toolbarPanel.add(toolbars[n]);
         }
         //        setMenubarState(toolbars[nToolbars-1], true);
         toolbarPanel.validate();
         toolbarPanel.repaint();
         validate();
         repaint();
     }

     public void addToolbar(StsToolbar toolbar_, int position)
     {
         final StsToolbar toolbar = toolbar_;
         if (hasToolbarNamed(toolbar.getName()))
         {
             StsToolkit.runLaterOnEventThread
             (
                 new Runnable()
                 {
                     public void run()
                     {
                        toolbar.setVisible(true);
                        setToolbarMenuState(toolbar, true);
                        repaintToolbar();
                     }
                 }
             );
             return;
         }
         if(Main.viewerOnly)
             toolbar.forViewOnly();
         toolbars = (StsToolbar[]) StsMath.arrayAddElement(toolbars, toolbar);
         StsToolkit.runLaterOnEventThread
         (
             new Runnable()
             {
                 public void run()
                 {
                    toolbarPanel.add(toolbar);
                    repaintToolbar();
                 }
             }
         );
     }

     public void closeToolbar(StsToolbar toolbar)
     {
         // toolbars = (StsToolbar[])StsMath.arrayDeleteElement(toolbars, toolbar);
         // toolbarPanel.remove(toolbar);
         toolbar.setVisible(false);
         //        toolbar.setVisible(false);
         setToolbarMenuState(toolbar, false);
         repaintToolbar();
     }

     private void repaintToolbar()
     {
         StsToolkit.runLaterOnEventThread
         (
             new Runnable()
             {
                 public void run()
                 {
                    toolbarPanel.validate();
                    toolbarPanel.repaint();
                    validate();
                    repaint();
                 }
             }
         );
     }

     private void setToolbarMenuState(StsToolbar toolbar, boolean state)
     {
         for (int n = 0; n < menubarToolbarItems.length; n++)
         {
             StsCheckboxMenuItem checkboxItem = menubarToolbarItems[n];
             if (!checkboxItem.getText().equals(toolbar.getName())) continue;
             checkboxItem.setSelected(state);
         }
     }

     public void removeToolbar(String toolbarName)
     {
         StsToolbar toolbar = StsToolbar.getToolbarNamed(toolbars, toolbarName);
         if (toolbar == null) return;
         closeToolbar(toolbar);
         repaint();
     }

     public StsViewGauges createGaugeView()
     {
         StsSensorClass sensorClass = (StsSensorClass)model.getStsClass("com.Sts.DBTypes.StsSensor");
         StsLiveWellClass liveWellClass = (StsLiveWellClass)model.getStsClass("com.Sts.DBTypes.StsLiveWell");
         StsObject[] objects =  new StsObject[0];
         if(sensorClass != null)
             objects = sensorClass.getSensors();
         StsObject[] wellObjects = new StsObject[0];
         if(liveWellClass != null)
             wellObjects = liveWellClass.getObjectList();

         if((objects.length == 0) && (wellObjects.length == 0))
         {
             new StsMessage(this, StsMessage.WARNING, "No gauge data was found in this project. Unable to launch gauge view.");
             return null;
         }

         StsViewGauges view = new StsViewGauges(this);
         if(view.getSelectedCurves() == null)
            return null;
         StsToolkit.centerComponentOnScreen(view);
         view.setVisible(true);

         // Add to window family
         getModel().addGaugeWindow(this, view);
         return view;
     }

     public StsViewTimeSeries createTimeSeriesView()
     {
         StsSensorClass sensorClass = (StsSensorClass)model.getStsClass("com.Sts.DBTypes.StsSensor");
         StsLiveWellClass liveWellClass = (StsLiveWellClass)model.getStsClass("com.Sts.DBTypes.StsLiveWell");
         StsObject[] objects =  new StsObject[0];
         if(sensorClass != null)
             objects = sensorClass.getSensors();
         StsObject[] wellObjects = new StsObject[0];
         if(liveWellClass != null)
             wellObjects = liveWellClass.getObjectList();

         if((objects.length == 0) && (wellObjects.length == 0))
         {
             new StsMessage(this, StsMessage.WARNING, "No time series data was found in this project. Unable to launch time series view.");
             return null;
         }
         StsViewTimeSeries view = null;
         if(!getModel().getProject().isRealtime())
            view = new StsViewTimeSeries(this);
         else
            view = new StsViewRealtimeSeries(this);

         if (view.selectedObjects == null || view.selectedObjects[0] == null)   // User canceled selection dialog
         {
             view = null;
             return null;
         }
         StsToolkit.centerComponentOnScreen(view);
         view.setVisible(true);

         // Add to window family
         getModel().addTimeWindow(this, view);

         // Reset the combo box to the value of the launching window.
         //        ((StsViewSelectToolbar)model.win3d.getToolbarNamed(StsViewSelectToolbar.NAME)).selectToolbarItem(currentView.name,true);
         StsViewSelectToolbar viewSelectToolbar = (StsViewSelectToolbar) getToolbarNamed(StsViewSelectToolbar.NAME);
         viewSelectToolbar.selectToolbarItem(StsViewTimeSeries.viewName, true);
         /*
         if (glPanel3d.getCurrentView() instanceof StsView3d)
             viewSelectToolbar.selectToolbarItem(StsView3d.viewClassname3d, true);
         else
             viewSelectToolbar.selectToolbarItem(StsViewCursor.viewClassname2d, true);
         //SAJ..............Need to reset selection to what it was not necessarily 3D view.
         */
         return view;
     }

     public void setView(StsViewItem viewItem)
     {

     }

     public void setView(String selectToolbarViewName)
     {
         StsViewSelectToolbar viewSelectToolbar = (StsViewSelectToolbar) getToolbarNamed(StsViewSelectToolbar.NAME);
         viewSelectToolbar.selectToolbarItem(selectToolbarViewName, true);
     }

     public synchronized void closeWindows()
     {
         model.closeWindows(this);
     }

     public void dispose()
     {
         if (currentViewPanel != null)
         {
             for (StsGLPanel3d glPanel3d : glPanel3ds)
                 glPanel3d.destroy();
         }
         super.dispose();
     }

     // Method: Help | About action performed

     public void helpAboutDialog()
     {
         StsWin3d_AboutBox dlg = new StsWin3d_AboutBox(this);
         Dimension dlgSize = dlg.getPreferredSize();
         Dimension frmSize = getSize();
         Point loc = getLocation();
         dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height - dlgSize.height) / 2 + loc.y);
         dlg.setModal(true);
         dlg.setVisible(true);
     }

     public void usageDialog()
     {
         StsWin3d_UsageBox dlg = new StsWin3d_UsageBox(this);
         Dimension dlgSize = dlg.getPreferredSize();
         Dimension frmSize = getSize();
         Point loc = getLocation();
         dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height - dlgSize.height) / 2 + loc.y);
         dlg.setModal(true);
         dlg.setVisible(true);
     }

     /*
         public void win3dDisplay()
         {
             if (glPanel3d != null)
             {
                 glPanel3d.repaint();
             }
         }
     */
     /** Base window doesn't have a 3d cursor, but this is needed for compatibility. */
     public void toggle3dCursor()
     {
     }

     public void adjustCursor(int dir, float dirCoor)
     {
         StsView[] displayedViews = getDisplayedViews();
         for (StsView view : displayedViews)
             view.adjustCursor(dir, dirCoor);
     }

     public void adjustCursorXY(float xCoor, float yCoor)
     {
         StsView[] displayedViews = getDisplayedViews();
         for (StsView view : displayedViews)
         {
             view.adjustCursor(StsCursor3d.XDIR, xCoor);
             view.adjustCursor(StsCursor3d.YDIR, yCoor);
         }
     }

     public void displayForeground()
     {
         if (StsLine.highlightedList != null)
             StsLine.highlightedList.display(getGlPanel3d());
     }

     /*
     public void refreshTabs()
     {
      model.refreshObjectPanel();
     }
     */
     /*
         public void repaint()
         {
             super.repaint();
         }
     */
     public boolean isMainWindow()
     {
         return this == getModel().win3d;
     }

     public JFrame getFrame()
     {
         return parent;
     }

     public StsWin3dBase[] getFamilyWindows()
     {
         return model.getWindowFamily(this).getWindows();
     }

     public StsWindowFamily getWindowFamily()
     {
         return model.getWindowFamily(this);
     }

     /*
     public void refreshTabs()
     {
      model.refreshObjectPanel();
     }
     */

     public void windowFamilyViewObjectChanged(Object source, Object object)
     {
         model.getWindowFamily(this).viewObjectChanged(source, object);
     }

     public void windowFamilyViewObjectRepaint(Object source, Object object)
     {
         model.getWindowFamily(this).viewObjectRepaint(source, object);
     }

     public void windowFamilyViewObjectChangedAndRepaint(Object source, Object object)
     {
         model.getWindowFamily(this).viewObjectChangedAndRepaint(source, object);
     }

     public void resetFamilyViews()
     {
         StsWin3dBase[] windows = getFamilyWindows();
         for (int i = 0; i < windows.length; i++)
         {
             if (windows[i].getGlPanel3d().getView() instanceof StsView2d)
             {
                 StsView2d view = (StsView2d) windows[i].getGlPanel3d().getView();
                 view.resetToOrigin();
                 break;
             }
         }
     }

     protected void addWindowCloseOperation()
     {
         setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
         addWindowListener
             (
                 new WindowAdapter()
                 {
                     public void windowClosing(WindowEvent e)
                     {
                         closeWindow();
                     }
                 }
             );
     }

     protected void closeWindow()
     {
         model.closeWindows(this);
     }

     protected void addMouseListener(StsWin3dBase window_)
     {
         final StsWin3dBase window = window_;
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

     static final long serialVersionUID = 1l;

     public void addSelectedObject(StsObject object)
     {
         selectedObjects = (StsObject[]) StsMath.arrayAddElementNoRepeat(selectedObjects, object);
     }

     public StsObject getSelectedObject(Class c)
     {
         if (selectedObjects == null) return null;
         for (int n = 0; n < selectedObjects.length; n++)
             if (selectedObjects[n].getClass() == c) return selectedObjects[n];
         return null;
     }

     public void saveGeometry()
     {
         if (currentViewPanel == null) return;
         size = getSize();
         location = getLocation();
         viewPanelSize = currentViewPanel.getSize();
         toolbarPanelSize = toolbarPanel.getSize();
         for (StsGLPanel3d glPanel3d : glPanel3ds)
             glPanel3d.saveGeometry();
     }

     public Point getLocation()
     {
         return new Point(getX(), getY());
     }

     public void componentResized(ComponentEvent e)
     {
         //       model.enableDisplay();
         //        StsException.systemDebug(this, "componentResized", "display.");
         super.validate();
         model.win3dDisplayAll(this);
     }

     /** Invoked when the component's position changes. */
     public void componentMoved(ComponentEvent e)
     {
         //       System.out.println("Frame moved.");
     }

     /** Invoked when the component has been made visible. */
     public void componentShown(ComponentEvent e)
     {

     }

     /** Invoked when the component has been made invisible. */
     public void componentHidden(ComponentEvent e)
     {

     }

     public void reshape(int x, int y, int width, int height)
     {
         //        System.out.println("Frame reshaped.");
         super.reshape(x, y, width, height);
     }

     public void displayPreviousObject(Class c)
     {
         StsObject prevObject = null;

         StsClass stsClass = model.getCreateStsClass(c);
         for (int i = 0; i < stsClass.getSize(); i++)
         {
             if (isMainWindow())
             {
                 prevObject = stsClass.setPreviousObject();
                 model.win3d.objectTreePanel.selected(prevObject);
             }
             else
             {
                 StsObject object = getCursor3d().getCurrentObject(c);
                 prevObject = stsClass.getPreviousObject(object);
             }
             if (prevObject != null && prevObject.getIsVisible())
                 break;
         }
         StsComboBoxToolbar toolbar = (StsComboBoxToolbar) model.win3d.getToolbarNamed(StsComboBoxToolbar.NAME);
         toolbar.comboBoxSetItem(prevObject);

         boolean changed = getCursor3d().setObject(prevObject);
         if (changed) repaint();
     }

     public void displayNextObject(Class c)
     {
         StsObject nextObject = null;

         StsClass stsClass = model.getCreateStsClass(c);
         for (int i = 0; i < stsClass.getSize(); i++)
         {
             if (isMainWindow())
             {
                 nextObject = stsClass.setNextObject();
                 model.win3d.objectTreePanel.selected(nextObject);
             }
             else
             {
                 StsObject object = getCursor3d().getCurrentObject(c);
                 nextObject = stsClass.getNextObject(object);
             }
             if (nextObject != null && nextObject instanceof StsMainObject)
             {
                 if (((StsMainObject) nextObject).getIsVisible()) break;
             }
         }
         if (nextObject == null) return;
         StsComboBoxToolbar toolbar = (StsComboBoxToolbar) getToolbarNamed(StsComboBoxToolbar.NAME);
         if (toolbar != null && toolbar.isVisible()) toolbar.comboBoxSetItem(nextObject);
         boolean changed = getCursor3d().setObject(nextObject);
         if (changed) repaint();
     }

     public void addGLPanel3d(StsGLPanel3d glPanel3d)
     {

     }

     public StsView getCurrentView() { return getGlPanel3d().view; }

     /*
         public StsGLPanel3d getSingleViewGlPanel3d()
         {
             if(currentViewPanel == null)
             {
                 currentViewPanel = singleViewPanel;
                 currentViewPanelIndex = -1;
             }
             else if(currentViewPanel != singleViewPanel)
             {
                 StsException.systemError(this, "getGlPanel3d", "Cannot get glPanel3d from multiViewPanel");
                 return null;
             }
             return singleViewPanel.glPanel3d;
         }
     */

     public StsView[] getDisplayedViews()
     {
         int nViews = 0;
         if (glPanel3ds != null) nViews = glPanel3ds.length;
         StsView[] views = new StsView[nViews];
         for (int n = 0; n < nViews; n++)
             views[n] = glPanel3ds[n].getView();
         return views;
     }

     public StsGLPanel3d[] getDisplayedGLPanels()
     {
         return glPanel3ds;
     }

     public void setViewsDefaultMouseModes()
     {
         if (currentViewPanel == null) return;
         StsView[] views = getDisplayedViews();
         for (StsView view : views)
             view.setDefaultMouseMode();
     }

     /**
      * If the multiViewPanel is in the window and this view is shown there, then we are ok...return.
      * If not, switch to singleViewPanel and change to or add this view.
      */

     public boolean hasView(Class viewClass)
     {
         for (StsGLPanel3d glPanel3d : glPanel3ds)
             if (glPanel3d.view.getClass() == viewClass) return true;
         return false;
     }

     public StsView getView(Class viewClass)
     {
         for (StsGLPanel3d glPanel3d : glPanel3ds)
             if (glPanel3d.view.getClass() == viewClass) return glPanel3d.view;
         return null;
     }

     public StsView[] getViews()
     {
         StsView[] views = new StsView[glPanel3ds.length];
         for(int n = 0; n < views.length; n++)
            views[n] = glPanel3ds[n].view;
         return views;
     }

     public StsGLPanel3d getGlPanelWithView(Class viewClass)
     {
         for (StsGLPanel3d glPanel3d : glPanel3ds)
             if (glPanel3d.view.getClass() == viewClass) return glPanel3d;
         return null;
     }
     /*
         public StsView getSinglePanelView(Class viewClass, String viewName)
         {
             return singleViewPanel.getView(viewClass, viewName);
         }

         public StsView setSinglePanelView(Class viewClass, String viewName)
         {
             return singleViewPanel.getView(viewClass, viewName);
         }
     */

     public StsView[] getViewsOfType(Class viewClass)
     {
         StsView[] viewsOfType = new StsView[0];
         StsView[] views = getDisplayedViews();
         for (StsView view : views)
             if (view.getClass() == viewClass)
                 viewsOfType = (StsView[]) StsMath.arrayAddElement(viewsOfType, view);
         return viewsOfType;
     }

     public void setDefaultViews()
     {
         Iterator<StsView> viewIterator = getViewIterator();
         while(viewIterator.hasNext())
         {
             StsView view = viewIterator.next();
             view.setDefaultView();
         }
     }

     public Iterator<StsView> getViewIterator()
     {
         return new ViewIterator();
     }

     private class ViewIterator implements Iterator<StsView>
     {
         int nViews = 0;
         int nView = 0;
         StsView[] views;
         StsView view;

         public ViewIterator()
         {
             views = getDisplayedViews();
             nViews = views.length;
         }

         public boolean hasNext()
         {
             if(nView >= nViews) return false;
             view = views[nView++];
             return true;
         }

         public StsView next()
         {
             return view;
         }

         public void remove()
         {
         }
     }

     public Iterator getViewIteratorOfType(Class viewType)
     {
         return new ViewIteratorOfType(viewType);

     }

     private class ViewIteratorOfType extends ViewIterator
     {
         Class viewType;

         public ViewIteratorOfType(Class viewType)
         {
             super();
             this.viewType = viewType;
         }

         public boolean hasNext()
         {
             while(nView < nViews)
             {
                 view = views[nView++];
                 if(viewType.isInstance(view)) return true;
             }
             return false;
         }
     }

     public void setDefaultView()
     {
         //        if(currentViewPanel == singleViewPanel)
         //            getGlPanel3d().setDefaultView();
     }

     public void initializeDefaultAction()
     {
         StsView[] views = getDisplayedViews();
         for (StsView view : views)
             view.initializeDefaultAction();
     }

     public void viewChanged()
     {
         StsView[] views = getDisplayedViews();
         for (StsView view : views)
             view.glPanel3d.viewChanged();
     }

     public void win3dDisplay()
     {
         StsView[] views = getDisplayedViews();
         for (StsView view : views)
             view.glPanel3d.repaint();
     }

     public boolean constructViewPanel(StsViewItem viewItem)
     {
         Class[] viewClasses = viewItem.viewClasses;
         int nViews = viewClasses.length;
         if (nViews == 0) return false;
         return constructViewPanel(viewItem.shortName, viewClasses);
     }

     public boolean checkConstructViewPanel(String viewName, Class viewClass)
     {
         if (hasView(viewClass)) return true;
         return constructViewPanel(viewName, viewClass);
     }

     public boolean checkConstructViewPanel(String viewName, StsView view)
     {
         if (hasView(view.getClass())) return true;
         return constructViewPanel(viewName, new StsView[] { view } );
     }
 
     public boolean constructViewPanel(String viewName, Class viewClass)
     {
         return constructViewPanel(viewName, new Class[]{viewClass});
     }

     public boolean constructViewPanel(String viewName, Class[] viewClasses)
     {
		 if(checkCurrentViewPanel(viewName)) return true;
         initializeGlPanel3ds(viewClasses);
		 return constructViewPanel(viewName);
	 }

     public boolean constructViewPanel(String viewName, StsView[] views)
     {
         if(checkCurrentViewPanel(viewName)) return true;
         initializeGlPanel3ds(views);
         return constructViewPanel(viewName);
     }

	 /** If currentViewPanel exists and has the same view, return true.
	  *  If it doesn't exist, return false.
	  *  If it exists but doesn't have the same view, remove it and return false.
	  */
	 private boolean checkCurrentViewPanel(String viewName)
	 {
         if (currentViewPanel == null) return false;
		 if(currentViewPanel.viewName == viewName) return true;
		 remove(currentViewPanel);
		 return false;
	 }

	 private boolean constructViewPanel(String viewName)
	 {
         StsMultiViewPanel currentViewPanel = StsMultiViewPanel.constructor(viewName, this);
		 if(currentViewPanel == null) return false;
		 this.currentViewPanel = currentViewPanel;
         return true;
     }
     /**
      * glPanel3ds is a set of gl panels which are reused.  Initialize the currentIf we currently have more than we need (nCurrentPanels > nPanels),
      * destroy the extra ones and trim glPanel3ds down to the ones being used.  If less that what we need, create new ones
      */
     private boolean initializeGlPanel3ds(Class[] viewClasses)
     {
         int nPanels = viewClasses.length;
         if (nPanels == 0)
         {
             StsException.systemError(this, "getCreateGlPanels", "Must have at least one panel for window.");
             return false;
         }
         constructGlPanel3ds(nPanels);
         // Initialize the view on each panel
         for (int n = 0; n < nPanels; n++)
            checkAddView(glPanel3ds[n], viewClasses[n]);
         return true;
     }
     
     private boolean initializeGlPanel3ds(StsView[] views)
     {
         int nPanels = views.length;
         if (nPanels == 0)
         {
             StsException.systemError(this, "getCreateGlPanels", "Must have at least one panel for window.");
             return false;
         }
         constructGlPanel3ds(nPanels);
         // Initialize the view on each panel
         for (int n = 0; n < nPanels; n++)
            checkAddView(glPanel3ds[n], views[n]);
         return true;
     }
     
     private void constructGlPanel3ds(int nPanels)
     {
         GraphicsDevice g = getGraphicsDevice();
         int nCurrentPanels = glPanel3ds.length;
         // if we have fewer number of panels than needed, reuse the ones we've got and add the new ones
         if (nCurrentPanels < nPanels)
         {
             // create an empty array of new panels
             StsGLPanel3d[] newGlPanel3ds = new StsGLPanel3d[nPanels];
             // copy the old panels to the new panels
             System.arraycopy(glPanel3ds, 0, newGlPanel3ds, 0, nCurrentPanels);
             // create the required number of added panels
             for (int n = nCurrentPanels; n < nPanels; n++)
                 newGlPanel3ds[n] = new StsGLPanel3d(model, actionManager, this);
             glPanel3ds = newGlPanel3ds;
         }
         // if we have more than the required number of panels, destroy the excess panels and reuse the remainder
         else if (nCurrentPanels > nPanels)
         {
             for (int n = nPanels; n < nCurrentPanels; n++)
                 glPanel3ds[n].destroy();
             glPanel3ds = (StsGLPanel3d[]) StsMath.trimArray(glPanel3ds, nPanels);
         }
     }

     private void checkAddView(StsGLPanel3d glPanel3d, Class viewClass)
     {
        glPanel3d.checkAddView(viewClass);
        glPanel3d.setSize(size00);
        glPanel3d.panelViewChanged = true;
     }
 
     private void checkAddView(StsGLPanel3d glPanel3d, StsView view)
     {
        glPanel3d.setView(view);
        glPanel3d.setSize(size00);
        glPanel3d.panelViewChanged = true;
     }
     
     public boolean initializeViewPanel()
     {
         int nPanels = glPanel3ds.length;
         if (nPanels == 0) return true;
         String[] viewNames = new String[nPanels];
         for (int n = 0; n < nPanels; n++)
         {
             glPanel3ds[n].initializeTransients(model, this);
             viewNames[n] = glPanel3ds[n].view.getViewName();
         }
         currentViewPanel.layoutPanel(viewNames);
         currentViewPanel.addComponents(glPanel3ds, viewNames);
         return true;
     }

     public boolean checkCursor3dIsOn() { return getCursor3d().isDisplay3dCursor(); }

     public StsCursor3d getCursor3d() { return getParentWindow().getCursor3d(); }

     public StsWin3dFull getParentWindow() { return model.getParentWindow(this); }

     public void checkAddDisplayableSections()
     {
         StsArrayList displayableClasses = model.getCursor3dTextureDisplayableClasses();
         for (int n = 0; n < displayableClasses.size(); n++)
         {
             Object classObject = displayableClasses.get(n);
             StsClass displayableClass = (StsClass) classObject;
             if (!displayableClass.hasObjects()) continue;
             StsClassCursor3dTextureDisplayable textureDisplayableClass = (StsClassCursor3dTextureDisplayable) classObject;
             for (int dir = 0; dir < 3; dir++)
                 checkAddDisplayableSection(textureDisplayableClass, dir);
         }
     }

     /**
      * array of displayableSections should correspond to array of displayableClasses; i.e., if displayableClasses(1) is StsSeismicVolumeClass,
      * then displayableSections(1) should be an instance of StsSeismicCursorSection.  If not, insert a new instance into the displayableSections
      * array and return it.
      */
     public void checkAddDisplayableSection(StsClassCursor3dTextureDisplayable displayableClass, int dir)
     {
         getCursor3d().checkAddDisplayableSection(displayableClass, dir);
     }

     public StsGLPanel3d getGlPanel3d()
     {
         if (glPanel3ds.length == 0) return null;
         return glPanel3ds[0];
     }

     public StsGLPanel3d[] getGLPanels()
     {
         return glPanel3ds;
     }

     public void toggleViewShift(ItemEvent e)
     {
         boolean viewShiftOn = (e.getStateChange() == ItemEvent.SELECTED);
         toggleViewShift(viewShiftOn);
     }

     public void toggleViewShift(boolean viewShiftOn)
     {
         this.viewShiftOn = viewShiftOn;
         for (StsGLPanel3d glPanel3d : getGLPanels())
             glPanel3d.changeModelView();
     }

     public void changeView(StsViewItem viewItem)
     {
         model.disableDisplay();
         constructViewPanel(viewItem);
         replaceViewPanel();
         validate();
         // validate();
         // invalidate();
         //        StsTextureList.clearAllTextures();
         model.enableDisplay();
         initAndStartGL();
     }

     protected void replaceViewPanel()
     {
         addBorderComponent(currentViewPanel, BorderLayout.CENTER);
     }
     /*
         public void startGL()
         {
             StsGLPanel3d[] glPanel3ds = this.currentViewPanel.getGLPanels();
             for(StsGLPanel3d glPanel3d : glPanel3ds)
                 glPanel3d.startGL();
         }
     */
 }
