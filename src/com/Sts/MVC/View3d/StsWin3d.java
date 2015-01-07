//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC.View3d;

import com.Sts.DB.*;
import com.Sts.Help.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;
import com.Sts.Workflow.*;
import com.Sts.WorkflowPlugIn.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.beans.*;
import java.lang.reflect.*;

public class StsWin3d extends StsWin3dFull implements StsSerializable
 {
     //    public boolean isCursor3dDisplayed = true;
     /** projection for all StsView3d views isVisible in windows; perspective if true; otherwise orthogonal */
     public boolean isPerspective = true;
     //    public StsCursor3d cursor3d = null;
     transient JSplitPane viewLogSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
     transient JSplitPane viewTreeSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);

     transient JTabbedPane tabbedPanels = new JTabbedPane();
     transient JTabbedPane tabbedTextPanels = new JTabbedPane();
     transient StsJPanel logPanel;

     transient public Sts3dMenuBar stsMenuBar;
     transient public StsObjectTree objectTreePanel;
     transient protected StsWorkflowPanel workflowPanel;

     // singleton: maintained thru window and model changes
     transient public StsSplashPanel splashPanel;

     transient public static final String VIEW_HELP_TARGET = "ui.3dDisplay.view";
     transient public static final String CURSOR_HELP_TARGET = "ui.3dDisplay.cursor";
     transient public static final String CONTROLS_HELP_TARGET = "ui.3dDisplay.controls";
     transient public static final String TREE_HELP_TARGET = "ui.tree.intro";
     transient public static final String STATUS_HELP_TARGET = "ui.status";

     public static final int LEFT = 0;
     public static final int RIGHT = 1;

     public StsWin3d()
     {
         // null constructor for serialization
     }

     public StsWin3d(StsModel model)
     {
         this.model = model;
         model.win3d = this;
         initializeActionManager();
         model.addWindowFamily(this);
         parentWindow = this;
         addComponentListener(this);
         windowName = "  Main Window";
     }

     static public StsWin3d constructDefaultMainWindow(StsModel model)
     {
         StsWin3d win3d = new StsWin3d(model);
         win3d.constructDefaultMainWindow();
         return win3d;
     }

     static public StsWin3d constructSplashMainWindow(StsModel model)
     {
         StsWin3d win3d = new StsWin3d(model);
         win3d.constructSplashWindow(model);
         return win3d;
     }

     public boolean start()
     {
         try
         {
             //			win3dInit();
             setEnabledMenuItems(true);
             model.win3d = this;
             model.setGlPanel3d(getGlPanel3d());
             return super.start();
         }
         catch (Exception e)
         {
             e.printStackTrace();
             return false;
         }
     }

     protected void constructSplashWindow(StsModel model)
     {
         constructBasicMainWindow();
         constructSplashPanel();
         initializeSplashLayout();
         start();
         // setSize(initialSize);
         //        centerComponentOnScreen();
     }

     protected void constructWindow()
     {
         constructBasicMainWindow();
         constructGraphicsPanel();
         initializeDefaultAction();
     }

     private void constructDefaultMainWindow()
     {
         constructBasicMainWindow();
         constructDefaultViewPanel();
         constructGraphicsPanel();
         initializeWindowLayout();
         initializeDefaultAction();
         start();
     }

     private void constructBasicMainWindow()
     {
         try
         {
             initializeDisplayIcon();
             setTitle("Version: " + StsModel.version + "  Project: " + model.getName() + windowName);

             addWindowCloseOperation(); // respond directly to window ops
             addMouseListener(this);  // used to catch window resize drag events
             // Container layout

             //initializeGridBagLayout();

             statusArea = new StsStatusArea(getActionManager()); // classInitialize the status area so listeners can connect

             HelpManager.deleteManager();
             stsMenuBar = new Sts3dMenuBar(model, this);
             setJMenuBar(stsMenuBar);

             objectTreePanel = new StsObjectTree(model, tabbedPanelSize);

             tabbedTextPanels.setTabPlacement(JTabbedPane.BOTTOM);

             //            statusArea.setPreferredSize(statusAreaSize);
             //            statusArea.setMinimumSize(new Dimension(200, 25));

             StsWorkflowPlugIn selectedPlugIn = model.getWorkflowPlugIn();
             if (selectedPlugIn == null)
             {
                 ClassLoader classLoader = getClass().getClassLoader();
                 Class plugInClass = classLoader.loadClass(model.workflowPlugInNames[0]);
                 Constructor constructor = plugInClass.getConstructor(new Class[0]);
                 selectedPlugIn = (StsWorkflowPlugIn) constructor.newInstance(new Class[0]);
             }
             workflowPanel = new StsWorkflowPanel(model, getActionManager(), selectedPlugIn, null);
             workflowPanel.setWorkflowPlugin(model.workflowPlugInNames[0]);

             tabbedPanels.add("Workflow", workflowPanel);
             tabbedPanels.add("Objects", objectTreePanel);

             //            pane.add(multiSplitPane, BorderLayout.CENTER);

             // add help links
             HelpManager.setContextSensitiveHelp(objectTreePanel, HelpManager.GENERAL, TREE_HELP_TARGET, this);
             HelpManager.setContextSensitiveHelp(getGlPanel3d(), HelpManager.GENERAL, VIEW_HELP_TARGET, this);
             HelpManager.setContextSensitiveHelp(statusArea, HelpManager.GENERAL, STATUS_HELP_TARGET, this);
             //			HelpManager.setContextSensitiveHelp(slider3dPanel, CONTROLS_HELP_TARGET, this);

             stsMenuBar.setEnabledMenuItems(false); // gray out everything but startup options
             //			setVisible(true);
         }
         catch (Exception e)
         {
             StsException.outputException("StsWin3d.classInitialize() failed.",
                 e, StsException.WARNING);
         }
     }

     public void initializeWindowLayout()
     {
//         initializeMultiViewSplitPane();
         viewTreeSplitPane.setOneTouchExpandable(true);
         viewTreeSplitPane.setContinuousLayout(true);
         if(!Main.viewerOnly)
             viewTreeSplitPane.add(tabbedPanels, JSplitPane.LEFT);
         else
            viewTreeSplitPane.remove(tabbedPanels);
         viewTreeSplitPane.setResizeWeight(0.0); // right component gets all the extra space
//         viewTreeSplitPane.setDividerLocation(0.30);
         viewTreeSplitPane.add(currentViewPanel, JSplitPane.RIGHT);

         // viewLogSplitPane.setBackground(Color.BLACK);
         // currentViewPanel.setBackground(Color.BLACK);

         viewLogSplitPane.setOneTouchExpandable(true);
         viewLogSplitPane.setContinuousLayout(true);
         viewLogSplitPane.add(viewTreeSplitPane, JSplitPane.TOP);
         viewLogSplitPane.add(logPanel, JSplitPane.BOTTOM);

         viewLogSplitPane.setResizeWeight(1.0); // top component gets all the extra space
         viewLogSplitPane.addPropertyChangeListener(new SplitChangeListener());

         // viewLogSplitPane.setBackground(Color.BLACK);

         initializeBorderLayout();
         addBorderComponent(toolbarPanel, BorderLayout.NORTH);
         addBorderComponent(viewLogSplitPane, BorderLayout.CENTER);

         tabbedPanels.setPreferredSize(tabbedPanelSize);
         tabbedTextPanels.setPreferredSize(tabbedTextPanelSize);
         tabbedPanels.setPreferredSize(tabbedPanelSize);
         // toolbarPanel.setPreferredSize(toolbarPanelSize);
         setPreferredSize(size);

         tabbedPanels.setMinimumSize(new Dimension(50, 0));
         currentViewPanel.setMinimumSize(size00);
         statusArea.setMinimumSize(new Dimension(50, 25));
         cursor3dPanel.setMinimumSize(new Dimension(50, 75));

//         viewTreeSplitPane.resetToPreferredSizes();
         viewLogSplitPane.resetToPreferredSizes();
     }

     private void constructSplashPanel()
     {
         splashPanel = new StsSplashPanel();
     }

     private void initializeSplashLayout()
     {
         viewTreeSplitPane.setOneTouchExpandable(true);
         viewTreeSplitPane.setContinuousLayout(true);
         viewTreeSplitPane.add(tabbedPanels, JSplitPane.LEFT);

         viewTreeSplitPane.setResizeWeight(1.0); // right component gets all the extra space
//         viewTreeSplitPane.setDividerLocation(0.30);
         viewTreeSplitPane.add(splashPanel, JSplitPane.RIGHT);

         initializeBorderLayout();
         addBorderComponent(viewTreeSplitPane, BorderLayout.CENTER);

         Dimension splashTabbedPanelSize = new Dimension(tabbedPanelSize.width, splashPanel.getHeight());
         tabbedPanels.setPreferredSize(splashTabbedPanelSize);
     }

     /**
      * A multisplitPane is used so we can independing adjust dividers between key window components.
      * The multiSplitPane is divided into  workflow and view panels in the top row and the log panel in the bottom row.
      * The log panel contains the 3dCursor panel, the tabbed message panels, and the statusArea panel.
      */
     private void constructGraphicsPanel()
     {
         splashPanel = null;
         constructTabbedMessagePanel();
         checkCreateToolbars();
         initializeTimeActionToolbar();
         constructCursor3d();
         constructLogPanel();
     }

     private void constructTabbedMessagePanel()
     {
         StsMessageFiles.initialize(model, tabbedTextPanels, model.getProject().getRootDirString());
         tabbedTextPanels.setSelectedIndex(0);
     }

     private void constructLogPanel()
     {
         logPanel = StsJPanel.noInsets();
         logPanel.gbc.fill = GridBagConstraints.BOTH;
         logPanel.gbc.weighty = 0.0;
         logPanel.add(cursor3dPanel);
         logPanel.gbc.weighty = 1.0;
         logPanel.add(tabbedTextPanels);
         logPanel.gbc.weighty = 0.0;
         logPanel.add(statusArea);
     }

     public void constructDefaultViewPanel()
     {
         if(!constructViewPanel(StsView3d.shortViewName3d, StsView3d.class)) return;
     }

     public void checkAddTimeActionToolbar()
     {
         if (hasToolbarNamed(StsTimeActionToolbar.NAME)) return;
         StsTimeActionToolbar timeActionToolbar = new StsTimeActionToolbar(this);
         addToolbar(timeActionToolbar);
         initializeTimeActionToolbar();
     }

     protected void createAddToolbars()
     {
         super.createAddToolbars();
         addToolbar(new StsMediaToolbar(this));
         addToolbar(new StsCollaborationToolbar(this));
     }

     protected void initializeTimeActionToolbar()
     {
         StsTimeActionToolbar toolbar = getTimeActionToolbar();
         if (toolbar == null) return;
         toolbar.initializeAction();
     }

     public void openWorkflowNewPage()
     {
         new StsMessage(this, StsMessage.INFO, "Creating a New Workflow currently unavailable.");
         //        URL editorPageURL = StsHtml.getHtml("workflowNew.html");
         //        if(workflowEditorPanel == null) workflowEditorPanel = StsWorkflowNew.constructor(editorPageURL, this);
         //        viewTreeSplitPane.remove(workflowSelectPanel);
         //        viewTreeSplitPane.add(workflowEditorPanel, JSplitPane.RIGHT);
     }

     public void openWorkflowEditPage()
     {
         new StsMessage(this, StsMessage.INFO, "Creating a New Workflow currently unavailable.");
         //        URL editorPageURL = StsHtml.getHtml("workflowEditor.html");
         //        if(workflowEditorPanel == null) workflowEditorPanel = StsWorkflowEditor.constructor(editorPageURL, this);
         //        viewTreeSplitPane.remove(workflowSelectPanel);
         //        viewTreeSplitPane.add(workflowEditorPanel, JSplitPane.RIGHT);
     }

     public void openWorkflowSavedPage()
     {
         new StsMessage(this, StsMessage.INFO, "Saving a Workflow currently unavailable.");
     }

     public void submitWorkflow()
     {
         //		viewTreeSplitPane.remove(workflowEditorPanel);
         //		viewTreeSplitPane.add(workflowSelectPanel, JSplitPane.RIGHT);
     }

     public void rebuildWorkflow(StsWorkflowPlugIn selectedPlugIn)
     {
         if ((selectedPlugIn == null) || (workflowPanel == null))
             return;
         tabbedPanels.remove(workflowPanel);
         tabbedPanels.remove(objectTreePanel);
         workflowPanel = null;
         workflowPanel = new StsWorkflowPanel(model, getActionManager(), selectedPlugIn, null);
         //        selectedPlugIn.addAdditionalToolbars(this);
         tabbedPanels.add("Workflow", workflowPanel);
         tabbedPanels.add("Objects", objectTreePanel);
         tabbedPanels.setSelectedIndex(0);
     }

     protected void createToolbarMenuItems()
     {
         super.createToolbarMenuItems();
         addToolbarMenuItems();
     }

     public void addToolbarMenuItems()
     {
         JMenu menuToolbars = stsMenuBar.menuToolbars;
         for (int n = 0; n < toolbars.length; n++)
             menuToolbars.add(menubarToolbarItems[n]);
     }

     public StsWorkflowPanel getWorkflowPanel()
     {
         return workflowPanel;
     }

     public void toggleVisibilityPane(ItemEvent e)
     {
         /*
         if(e.getStateChange() == ItemEvent.SELECTED)
         {
             if(lastDividerLocation != -1)
                 viewTreeSplitPane.setDividerLocation(lastDividerLocation);
             else
                 viewTreeSplitPane.setDividerLocation(0.30);
         }
         else
         {
             lastDividerLocation = viewTreeSplitPane.getDividerLocation();
             viewTreeSplitPane.setDividerLocation(0.0);
         }
         */
     }

     public StsTimeActionToolbar getTimeActionToolbar()
     {
         return (StsTimeActionToolbar) getToolbarNamed(StsTimeActionToolbar.NAME);
     }

     // Method: win3dInit

     private void win3dInit() throws Exception
     {
         initializeProperties();
         //        initializeToolbars();
     }

     public void resetCursorPanel()
     {
         StsCursor3d cursor3d = getCursor3d();
         cursor3d.initialize();
         //        if(cursor3dPanel != null)  this.cursor3dPanel.setValues();
     }

     /** Initialize any properties, flags after win3d is initialized */
     private void initializeProperties()
     {
         boolean useDisplayLists = model.useDisplayLists;
         stsMenuBar.setCheckboxUseDisplayLists(useDisplayLists);
     }

     public void setEnabledMenuItems(boolean enabled)
     {
         if (stsMenuBar != null) stsMenuBar.setEnabledMenuItems(enabled);
     }

     public StsMenuItem getWorkflowHelpButton()
     {
         return stsMenuBar.menuStepsManual;
     }
     public void enable3dCursor(boolean enabled)
     {
         if (isCursor3dDisplayed == enabled) return;
         isCursor3dDisplayed = enabled;
         this.cursor3dPanel.setEditable(enabled);
         setDisplay3dCursorPanel(isCursor3dDisplayed);
     }

     public void setDisplay3dCursorPanel(boolean selected)
     {
         HelpManager.setContextSensitiveHelp(cursor3dPanel, HelpManager.GENERAL, CURSOR_HELP_TARGET, this);
         cursor3dPanel.setVisible(selected);
         setDisplay3dCursor(selected);
     }

     public StsCursor3dPanel getCursor3dPanel()
     {
         return cursor3dPanel;
     }

     public void setZDomain()
     {
         if (model != null)
         {
             if (model.getProject() != null)
             {
                 model.getProject().setZDomain();
                 //               defaultToolbar.changeTimeDepthIcon(model.getProject().getZDomain());
             }
         }
     }

     public void repaint()
     {
         if (currentViewPanel != null) currentViewPanel.repaint();
         setDividerState();
         tabbedPanels.validate();
     }

     private void setDividerState()
     {
         /*
         int maxLoc = viewTreeSplitPane.getMaximumDividerLocation();
         boolean selected = viewTreeSplitPane.getDividerLocation() <= maxLoc;
 //		maxLoc -= (viewTreeSplitPane.getDividerSize() + 2);
 //        System.out.println("Divider changed. maxLoc: " + maxLoc +
 //                            " current: " + viewTreeSplitPane.getDividerLocation());
 //		boolean selected = false;
         if(stsMenuBar == null) return;
         ((Sts3dMenuBar)stsMenuBar).menuVisibilityDetail.setState(selected);
         */
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

     public void rebuildMenubar()
     {
         stsMenuBar.buildMenuDatabases();
     }

     protected void closeWindow()
     {
         model.appClose();
     }

     static final long serialVersionUID = 1l;
/*
     public void toggleProjection(ItemEvent e)
     {
         boolean isPerspective = (e.getStateChange() == ItemEvent.SELECTED);
         toggleProjection(isPerspective);
     }

     public void toggleProjection(boolean isPerspective)
     {
         if(model.getProject().getIsPerspective() == isPerspective) return;
         model.getProject().setIsPerspective(isPerspective);
         toggleProjection();
     }
  */
     public void toggleProjection()
     {
         StsWin3dBase[] windowsWithView = model.getWindowsWithViewOfType(StsView3d.class);
         for (int n = 0; n < windowsWithView.length; n++)
         {
             StsView[] viewsOfType = windowsWithView[n].getViewsOfType(StsView3d.class);
             if (viewsOfType.length == 0) continue;
             for (StsView view : viewsOfType)
             {
                 StsGLPanel3d glPanel3d = ((StsView3d) view).glPanel3d;
                 glPanel3d.projectionChanged = true;
                 glPanel3d.repaint();
             }
         }
     }

     public void saveGeometry()
     {
         size = getSize();
         location = getLocation();
         viewPanelSize = currentViewPanel.getSize();
         tabbedPanelSize = tabbedPanels.getSize();
         for(StsGLPanel3d glPanel3d : glPanel3ds)
            glPanel3d.saveGeometry();
     }

    protected void replaceViewPanel()
    {
        viewTreeSplitPane.setRightComponent(currentViewPanel);
    }
 }
