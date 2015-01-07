//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC.View3d;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/** Independent window with 3d cursor, but no object panel. */
 public class StsWin3dFull extends StsWin3dBase implements StsSerializable
 {
     public StsCursor3d cursor3d = null;
     public boolean isCursor3dDisplayed = true;
     private boolean cursorDisplayXY = true;
     /** Panel displaying cursor controls: doesn't exist in StsWin3dBase (an auxiliary window); only in StsWin3d or StsWin3dFull */
     transient public StsCursor3dPanel cursor3dPanel = null;

     public StsWin3dFull()
     {
     }

     /** Main window for a new family of windows.  Spawned from original 3d window when toolbar button selected. */
     public StsWin3dFull(StsModel model)
     {
         this.model = model;
         actionManager = new StsActionManager(model);
         parentWindow = this;
         if (model.addWindowFamily(this) == false)
             closeWindow();

         try
         {
             windowName = model.getWindowFamilyTitle(this);
             constructWindow();
         }
         catch (Exception e)
         {
             e.printStackTrace();
         }
     }

     protected void constructWindow()
     {
         try
         {
             initializeDisplayIcon();
             setTitle(model.getName() + " " + windowName);
             addWindowCloseOperation(); // respond directly to window ops
             addMouseListener(this);  // used to catch window resize drag events
             checkCreateToolbars();
             constructCursor3d();
             //             centerComponentOnScreen();
         }
         catch (Exception e)
         {
             StsException.outputException("StsWin3d.classInitialize() failed.", e, StsException.WARNING);
         }
     }

     public void initializeActionManager()
     {
         actionManager = new StsActionManager(model);
     }

     public void initializeWindowLayout()
     {
         addBorderComponent(toolbarPanel, BorderLayout.NORTH);
         addBorderComponent(currentViewPanel, BorderLayout.CENTER);
         addBorderComponent(cursor3dPanel, BorderLayout.SOUTH);
         // currentViewPanel.setMinimumSize(new Dimension(0,0));
         currentViewPanel.setPreferredSize(viewPanelSize);
     }

     protected void constructCursor3d()
     {
         if (cursor3d == null)
         {
             cursor3d = new StsCursor3d(model, this);
             cursor3d.initialize();
         }
         else
         {
             cursor3d.model = model;
             cursor3d.window = this;
         }

         cursor3dPanel = new StsCursor3dPanel(model, this, cursor3d, cursor3dPanelSize);
         setDisplay3dCursorPanel(isCursor3dDisplayed);
     }

     protected void initializeCursor3dTransients(StsModel model, StsWin3dBase win3d)
     {
         cursor3d.initializeTransients(model, this);
     }

     protected void createAddToolbars()
     {
         super.createAddToolbars();
         if(!Main.viewerOnly)
         {
             addToolbar(new StsIntraFamilyToolbar(this));
             StsPreStackLineSet.checkAddToolbar(model, this);
         }
     }

     // Accessors

     //	public StsWindowActionManager getActionManager() { return actionManager; }

     public StsModel getModel()
     {
         return model;
     }

     /** Initialize any properties, flags after win3d is initialized */
     private void initializeProperties()
     {
     }

     public void toggle3dCursor()
     {
         isCursor3dDisplayed = !isCursor3dDisplayed;
         setDisplay3dCursorPanel(isCursor3dDisplayed);
     }

     public void setDisplay3dCursorPanel(boolean selected)
     {
         Container pane = this.getContentPane();
         if (selected)
             pane.add(cursor3dPanel, BorderLayout.SOUTH);
         else
             pane.remove(cursor3dPanel);
         pane.validate();

         setDisplay3dCursor(selected);
         win3dDisplay();
     }

     public void cursorPickSetup()
     {
         if (checkCursor3dIsOn() == false)
         {
             setDisplay3dCursorPanel(true);
             setDisplay3dCursor(true);
             //            win3dDisplay();
         }
     }

     public void repaint()
     {
         super.repaint();
     }

     protected void closeWindow()
     {
         model.closeFamily(this);
     }

     public void windowClosing(WindowEvent e)
     {
         model.closeFamily(this);
     }

     // Stubbed in to satisfy dependency until I figure out where they should be.
     public void trackCursor()
     {}

     public void trackMouse()
     {}

     public void trackUser()
     {}

     static final long serialVersionUID = 1l;

     /** If main window is iconified, icon all children as well */
     public void windowIconified(WindowEvent e)
     {
         model.iconifyWindowFamily(this, JFrame.ICONIFIED);
     }

     /** If main window is de-iconified, de-icon all children as well */
     public void windowDeiconified(WindowEvent e)
     {
         model.iconifyWindowFamily(this, JFrame.NORMAL);
     }

     public StsCursor3d getCursor3d() { return cursor3d; }

     public void adjustCursorAndSlider(int dir, float dirCoor)
     {
         adjustSlider(dir, dirCoor);
         adjustCursor(dir, dirCoor);
     }

     public void adjustSlider(int dir, float dirCoor)
     {
         if (cursor3dPanel != null) cursor3dPanel.setSliderValue(dir, dirCoor);
     }

     public void incrementCursor(int dir)
     {
         cursor3dPanel.step(dir, true);
     }
     public void decrementCursor(int dir)
     {
         cursor3dPanel.step(dir, false);
     }

     public boolean getCursorDisplayXY()
     {
         return cursorDisplayXY;
     }

     public void setCursorDisplayXY(boolean isXY)
     {
         this.cursorDisplayXY = isXY;
     }

     public void setCursorDisplayXYAndGridCheckbox(boolean isXY)
     {
         setCursorDisplayXY(isXY);
         if (cursor3dPanel != null) cursor3dPanel.setGridCheckboxState(isXY);
     }

     public void setDisplay3dCursor(boolean display3dCursor) { getCursor3d().setDisplay3dCursor(display3dCursor); }


     public void setZDomainSupported(byte zDomainSupported)
     {
         StsComboBoxFieldBean timeDepthBean = (StsComboBoxFieldBean) getToolbarComponentNamed(
             StsDefaultToolbar.NAME, StsDefaultToolbar.TIME_DEPTH);
         if (timeDepthBean == null) return;
         timeDepthBean.setEditable(zDomainSupported == StsProject.TD_TIME_DEPTH);
     }

     public void setZDomain(byte zDomain)
     {
         StsComboBoxFieldBean timeDepthBean = (StsComboBoxFieldBean) getToolbarComponentNamed(
             StsDefaultToolbar.NAME, StsDefaultToolbar.TIME_DEPTH);
         if (timeDepthBean == null) return;
         timeDepthBean.setSelectedItem(StsParameters.TD_ALL_STRINGS[zDomain]);
     }
 }
