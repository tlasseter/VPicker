package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.StsView;
import com.Sts.MVC.ViewWell.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import javax.media.opengl.glu.GLU;

/**
 * Sts Well Viewer
 * Description:  Well Model-Viewer
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 *
 * @author T.J.Lasseter
 * @version 1.0
 * @description wellViewModel includes the depth range of the well, list of logTracks being isVisible, and parameters required to build the view in StsWellWindowFrame
 */
abstract public class StsWellViewModel extends StsSerialize
{
    /** well this is a viewModel for */
    public StsWell well;
    /** well measured depth min */
    public double mdepthMin;
    /** well measured depth max */
    public double mdepthMax;
    /** well window is currently being isVisible */
    public boolean isVisible = true;
    /** set of isVisible panels objects (log, seismic, and VSP panels) for this well model */
    public StsWellView[] wellViews;
    /** current zoomLevel index */
    public int zoomCurrentLevel = 0;
    /** height of wellTrack, logCurve and vsp displayPanels; includes any boxes at the top of each wellView */
    public int displayHeight = defaultHeight;
    /** measured depth at bottom of window; top is computed on demand based on window size */
    private double windowMdepthMax;
    /** curve track width for this well model */
    public int curveTrackWidth;
    /** horizontal cursor pixel Y position (GL: up from bottom!) */
    public int cursorY;
    /** cursor Z position * */
    public double cursorMdepth;

    /** Action manager for this window, all panels */
    transient public StsActionManager actionManager;
    /** panel for the view of this well model */
    transient public StsWellWindowPanel wellWindowPanel;
    /** current log track */
    transient StsLogCurvesView currentCurvesView = null;
    /** track units: mdepth, depth, or time */
    transient String units = MDEPTH;
    /** current view that has focus */
    transient public StsWellView currentWellView = null;
    /** current zoomLevel */
    transient public StsWellWindowZoomLevel zoomLevel;
    /** cursor is currently being moved */
    transient public boolean cursorPicked = false;

    /** minimum zoomLevel index allowed (see StsWellWindowZoomLevel */
    transient protected int zoomMinLevel = 0;

    /** maximum zoomLevel index allowed */
    transient protected int zoomMaxLevel = 0;

    /** view of the well track: vertical line against depth scale from wellMinZ to wellMaxZ */
    transient public StsWellTrackView wellTrackView;

    /** current Y position of mouse measured up from bottom */
    transient protected int glMouseY;

    /** current Z value of mouse */
    transient protected double mouseZ;

    transient private int logsPerTrack = 2;

    transient int nSubWindows = 0; //

    /** array of zoomLevels defining pixels/tick, ticks/per label, and units/pixel, etc */
    static StsWellWindowZoomLevel[] zoomLevels = StsWellWindowZoomLevel.zoomLevels;

    /** number of possible zoomLevels */
    static int nZoomLevels = StsWellWindowZoomLevel.nZoomLevels;

    /** default window height */
    static protected int defaultHeight = 500;

    /** default panelWidth applied when new panel is constructed */
    static protected int defaultPanelWidth = 150;

    /** initial window X */
    static int windowX = 100;

    /** initial window Y */
    static int windowY = 100;

    /** window dX to next window */
    static int windowDX = 25;

    /** window dY to next window */
    static int windowDY = 25;

    /** unit types */
    static public String MDEPTH = "MDepth";
    static public String DEPTH = "Depth";
    static public String TIME = "Time";
    static public String MDEPTH_TIME = "MDepth & Time";
    static public String DEPTH_TIME = "Depth & Time";

    /** operation types */
    static public String NO_ACTION = "No action";
    static public String ADD_LOG_CURVES = "Add Log Panel";
    static public String ADD_SYNTHETIC = "Add Synthetic Seismic";
    static public String EDIT_LOG_CURVE = "Edit Log Curves";
    static public String ADD_WELL_TRACK = "Add Well Track";
    static public String ADD_TRACE_DISPLAY = "Add Trace Display";
    static public String ADD_TD_CURVE = "Add T/D Curve";
    static public String ADD_INTERVAL_VELOCITY = "Add Interval Velocity";
    static public String ADD_VSP = "Add VSP Panel";
    static public String DISPLAY_VSP = "Display VSP";
    static public String ADD_CURTAIN = "Add Seismic Curtain";
    static public String ADD_DTS_PANEL = "Add DTS Panel";
    static public String[] operationTypes = new String[]
        {
            NO_ACTION, ADD_LOG_CURVES, ADD_TRACE_DISPLAY, ADD_SYNTHETIC, ADD_CURTAIN, EDIT_LOG_CURVE, ADD_WELL_TRACK, ADD_TD_CURVE, ADD_INTERVAL_VELOCITY
        };
    static public String[] operationTypesVSP = new String[]
        {
            NO_ACTION, ADD_LOG_CURVES, /*EDIT_LOG_CURVE,  ADD_WELL_TRACK,*/ ADD_CURTAIN, ADD_VSP
        };

    /** spectrum used for well line colors */
//    static StsSpectrum colorSpectrum;
    /*
      static
      {
      try
       {
        colorSpectrum = new StsSpectrum("Basic");
        colorSpectrum.setBasic32Colors();
       }
       catch(Exception e)
       {
        System.out.println("StsWellViewModel static initializer for colorSpectrum failed.");
       }
      }
      */
    static public NumberFormat labelFormatter = NumberFormat.getNumberInstance();
    static public GLBitmapFont font;
    static public int fontHeight;
    static public int fontWidth;

    abstract public Frame getParentFrame();

    abstract protected Component getCenterComponent();

    abstract public void rebuild();

    public StsWellViewModel()
    {
    }

    public StsWellViewModel(StsWell well, StsModel model)
    {
        this.well = well;
        actionManager = new StsActionManager(model);
        layoutWellWindow();
    }

    public void initializeWellWindowPanel()
    {
        wellWindowPanel = new StsWellWindowPanel();

        if (wellViews != null) nSubWindows = wellViews.length;
        if (nSubWindows == 0)
        {
            wellTrackView = new StsWellTrackView();
            addWellView(wellTrackView);
        }
        for (int n = 0; n < nSubWindows; n++)
            wellViews[n].initializeView(this, currentModel, actionManager, n);
    }

    private void addWellView(StsWellView wellView)
    {
        wellViews = (StsWellView[]) StsMath.arrayAddElement(wellViews, wellView, StsWellView.class);
        nSubWindows++;

//        rebuild();
    }

    public boolean initialize(StsModel model)
    {
//		wellViewModelClass = (StsWellViewModelClass)currentModel.getCreateStsClass(StsWellViewModel.class);
        return true;
    }

    // jbw new method to construct font (requires active GL context)
    public static void initializeFont(GL gl)
    {
//		if(font == null)
        {
            font = GLHelvetica12BitmapFont.getInstance(gl);
            fontHeight = font.getCharacterHeight('A');
            fontWidth = font.getCharacterWidth('Z');
        }
    }

    public void addKeyListener(KeyListener listener)
    {
        wellWindowPanel.addKeyListener(listener);
//        if(wellWindowFrame != null) wellWindowFrame.addKeyListener(listener);
//        else if(wellWindowPanel != null) wellWindowPanel.addKeyListener(listener);
    }


    public String getName()
    {
        if (well == null) return new String("null");
        return well.getName();
    }

    public StsWell getWell()
    {
        return well;
    }

    public StsWellWindowPanel getWellWindowPanel()
    {
        return wellWindowPanel;
    }

    public void addToWellWindowPanel(StsGLPanel glPanel, Object constraints)
    {
        wellWindowPanel.innerPanel().add(glPanel, constraints);
        //GraphicsDevice gDevice = glPanel.getGraphicsConfiguration().getDevice();
        glPanel.initAndStartGL(glPanel.getGraphicsDevice());
        // wellWindowFrame.pack();
        // wellWindowFrame.repaint();
    }

    public void setHeight(int height)
    {
        displayHeight = height;
    }

    public void setCurveTrackSize(Dimension size)
    {
        /*
           curveTrackWidth = size.width;
           currentCurveTrackWidth = size.width;
           displayHeight = size.height;
           for (int n = 0; n < logTracks.getSize(); n++)
           {
            StsLogTrack logTrack = (StsLogTrack) logTracks.getElement(n);
            logTrack.getLogTrackView().setDisplayPanelSize(size);
           }
           */
    }

    public int getCurveTrackWidth()
    {
        return curveTrackWidth;
    }

    public int getDisplayHeight()
    {
        return displayHeight;
    }

    public void editLogCurves()
    {
        // Display log curve edit dialog
        StsWellEditDialog wellDialog = new StsWellEditDialog(this, well, false);
        wellDialog.setVisible(true);
    }

    public void runSelectLogCurves()
    {

        //Runnable runSelectLogCurves = new Runnable()
        //{
        //	public void run()
        //	{
        selectLogCurves();
        //	}
        //};

        //StsToolkit.runWaitOnEventThread(runSelectLogCurves);
        //Thread selectLogCurvesThread = new Thread(runSelectLogCurves);
        //selectLogCurvesThread.start();
//        runPickSurface.run();
    }

    public void selectLogCurves()
    {
        selectLogCurves(null);
    }

    public void selectLogCurves(StsLogCurvesView curvesView)
    {
        currentCurvesView = curvesView;

        StsLogCurve[] logCurves = getAvailableLogCurves();
        if (logCurves == null)
        {
            message(StsMessage.INFO, "No well curves available.");
            addSelectedLogCurves(new StsLogCurve[0]);
            return;
        }
        int nLogCurves = logCurves.length;
        if (nLogCurves == 0)
        {
            message(StsMessage.INFO, "No well curves available.");
            addSelectedLogCurves(new StsLogCurve[0]);
            return;
        }
        String[] logCurveNames = new String[nLogCurves];
        for (int n = 0; n < nLogCurves; n++)
        {
            logCurveNames[n] = logCurves[n].getName();
        }
        Object[] selectedLogCurves = StsListSelectionDialog.getMultiSelectFromListDialog(getParentFrame(), getCenterComponent(), "Log Curves",
            "Select log curves from list.", logCurveNames, logCurves);
        if (selectedLogCurves == null)
        {
            message(StsMessage.INFO, "No logCurve files selected.");
            return;
        }
        addSelectedLogCurves(selectedLogCurves);
    }

    public String selectLogCurvesForSynth(String type)
    {

        StsLogCurve[] logCurves = (StsLogCurve[]) well.getLogCurves().getCastListCopy();
        if (logCurves == null)
        {
            message(StsMessage.INFO, "No well curves available.");

            return null;
        }
        int nLogCurves = logCurves.length;
        if (nLogCurves == 0)
        {
            message(StsMessage.INFO, "No well curves available.");
            return null;
        }
        String[] logCurveNames = new String[nLogCurves];
        for (int n = 0; n < nLogCurves; n++)
        {
            logCurveNames[n] = logCurves[n].getName();
        }
        Object[] selectedLogCurves = StsListSelectionDialog.getMultiSelectFromListDialog(getParentFrame(), getCenterComponent(), "Log Curves",
            "Select " + type + " log curve from list for synthetic calculation.", logCurveNames, logCurves);
        if (selectedLogCurves == null)
        {
            message(StsMessage.INFO, "No logCurve files selected.");
            return null;
        }
        if (selectedLogCurves.length != 1)
        {
            message(StsMessage.INFO, "I need 1 log curve.");
            return null;
        }

        String list = new String();
        list = ((StsLogCurve) selectedLogCurves[0]).getName();

        return list;
    }

    private void addSelectedLogCurves(Object[] selectedLogCurves)
    {
        int nSelectedCurves = selectedLogCurves.length;

        StsLogCurve[] selectedCurves = new StsLogCurve[nSelectedCurves];
        for (int n = 0; n < nSelectedCurves; n++)
        {
            selectedCurves[n] = (StsLogCurve) selectedLogCurves[n];

        }
        addLogCurves(selectedCurves);
    }

    public void removeLogPanel(StsLogCurvesView curvesView)
    {
        curvesView.getLogTrack().removeAllLogCurves();
        curvesView.removePanel();
    }

    public void selectVSPs()
    {
        StsVsp[] vsps = getAvailableVSPs();
        int nv = vsps.length;
        if (nv == 0)
        {
            message(StsMessage.INFO, "No VSPs available or associated with this well.");
            return;
        }
        String[] names = new String[nv];
        for (int n = 0; n < nv; n++)
        {
            names[n] = vsps[n].getName();

        }
        Object[] selected = StsListSelectionDialog.getMultiSelectFromListDialog(getParentFrame(), getCenterComponent(), "VSPs",
            "Select VSPs from list.", names, vsps);
        if (selected == null)
        {
            message(StsMessage.INFO, "No VSPs selected.");
            return;
        }

        int nSelected = selected.length;

        StsVsp[] selectedVSP = new StsVsp[nSelected];
        for (int n = 0; n < nSelected; n++)
        {
            selectedVSP[n] = (StsVsp) selected[n];

        }
        addVSPs(selectedVSP);
    }

    private void message(int type, String message)
    {
//        if(wellWindowFrame != null)
        new StsMessage(getParentFrame(), type, message);
//        else if(wellWindowPanel != null)
//            new StsMessage(model.win3d, type, message);
    }

    public void doAddLogCurves(StsLogCurve[] logCurves)
    {
        if (logCurves == null) return;
        int nLogCurves = logCurves.length;
        if (nLogCurves == 0)
        {
            addLogCurve(null);
            return;
        }
        for (int n = 0; n < nLogCurves; n++)
            addLogCurve(logCurves[n]);
    }

    // all of this backhandedness is to overcome the parenting problem on multiple windows.
    // we need to get the adding of curves out of te SWING thread so that we can use SWING
    // events to manage adding to screen 1
    // jogl b5 fixes this issue - jbw
    public void addLogCurves(StsLogCurve[] logCurves)
    {
        doAddLogCurves(logCurves);
        //AddLogCurvesWorker t = new AddLogCurvesWorker(logCurves);
        //t.start();
    }

    public void removeLogCurves(StsLogCurve[] logCurves)
    {
        if (logCurves == null)
        {
            return;
        }
        int nLogCurves = logCurves.length;
        for (int n = 0; n < nLogCurves; n++)
            removeLogCurve(logCurves[n]);
    }

    public StsLogCurve[] getAvailableLogCurves()
    {
        StsLogCurve[] logCurves = (StsLogCurve[]) well.getLogCurves().getCastListCopy();
        if (logCurves == null || logCurves.length == 0) return logCurves;
        StsLogCurve[] displayedLogCurves = getDisplayedLogCurves();

        if (displayedLogCurves == null) return logCurves;
        int nDisplayedLogCurves = displayedLogCurves.length;
        if (nDisplayedLogCurves == 0) return logCurves;

        for (int d = 0; d < nDisplayedLogCurves; d++)
        {
            StsLogCurve displayedLogCurve = displayedLogCurves[d];
            logCurves = (StsLogCurve[]) StsMath.arrayDeleteElement(logCurves, displayedLogCurve);
        }
        return logCurves;
    }

    private StsLogCurve[] getDisplayedLogCurves()
    {
        int nViews = wellViews.length;
        StsLogCurve[] logCurves = new StsLogCurve[0];
        for (int n = 0; n < nViews; n++)
        {
            if (wellViews[n] instanceof StsLogCurvesView)
            {
                StsLogTrack logTrack = (StsLogTrack) wellViews[n].displayedObject;
                logCurves = (StsLogCurve[]) StsMath.arrayAddArray(logCurves, logTrack.logCurves.getElements());
            }
        }
        return logCurves;
    }

    public boolean isDisplayingLogCurve(StsLogCurve logCurve)
    {
        int nViews = wellViews.length;
        for (int n = 0; n < nViews; n++)
        {
            if (wellViews[n] instanceof StsLogCurvesView)
            {
                StsLogTrack logTrack = (StsLogTrack) wellViews[n].displayedObject;
                if(logTrack.isDisplayingLogCurve(logCurve))
                    return true;
            }
        }
        return false;
    }

    public StsLogCurvesView addLogCurve(StsLogCurve logCurve)
    {
        if(logCurve == null) return null;
        logCurve.checkLoadVectors();
        StsLogCurvesView logCurvesView = addCurveToLogTrack(logCurve);
        //StsToolkit.runWaitOnEventThread(new Runnable() { public void run() { StsWellViewModel.rebuild(); } } );
        rebuild();
        repaint();
        return logCurvesView;
    }

    public void removeLogCurve(StsLogCurve logCurve)
    {
        if(logCurve == null) return;
        removeCurveFromLogTrack(logCurve);
        rebuild();
        repaint();
    }

    /*
      public void addWellTrack()
      {
       JPanel panel = new JPanel();
       panel.setSize(new Dimension(125, 500));
       panel.setBackground(Color.BLUE);
       wellWindowPanel.add(panel, new GridBagConstraints(1 + nPanelsAdded, 1, 1, 1, 1.0, 1.0,
        GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
       nPanelsAdded++;
       rebuild();
 //        int nTracks = 1;
 //        StsWellTrackView wellTrackView = new StsWellTrackView(this, model, wellWindowPanel, 1);
      }
      */


    private StsLogTrack getLastLogTrack()
    {
        int nWellViews = wellViews.length;
        for (int n = nWellViews - 1; n > 0; n--)
        {
            if (wellViews[n] instanceof StsLogCurvesView)
                return (StsLogTrack) wellViews[n].displayedObject;
        }
        return null;
    }

    public StsLogCurvesView addCurveToLogTrack(StsLogCurve logCurve)
    {

        //int nTracks = displayedPanels.getSize();

        if (currentCurvesView == null || currentCurvesView.getNLogCurves() >= StsLogTrack.maxCurvesPerTrack)
        {
            currentCurvesView = addLogTrackToWindow();
//			nSubWindows++;

        }
        currentCurvesView.addLogCurve(logCurve);
        return currentCurvesView;
    }

    public void removeCurveFromLogTrack(StsLogCurve logCurve)
    {

        if (currentCurvesView == null)
            return;

        currentCurvesView.removeLogCurve(logCurve);
        return;
    }

    public StsLogCurvesView addLogTrackToWindow()
    {
        if (wellWindowPanel == null) return null;
        int width = wellWindowPanel.getWidth();
        int newWidth = width + curveTrackWidth;
        wellWindowPanel.setSize(newWidth, displayHeight);
        StsLogCurvesView view = new StsLogCurvesView(this, currentModel, actionManager, nSubWindows);
        addWellView(view);
        return view;
    }

    public void layoutWellWindow()
    {
        curveTrackWidth = defaultPanelWidth;
        //displayHeight = defaultDisplayHeight;
        mdepthMin = well.getMinMDepth();
        mdepthMax = well.getMaxMDepth();

        double wellDZ = mdepthMax - mdepthMin;
        double desiredPixelsPerUnit = (double) displayHeight / wellDZ;
        for (int n = nZoomLevels - 1; n > 0; n--)
        {
            if (zoomLevels[n].pixelsPerUnit < desiredPixelsPerUnit)
            {
                zoomMinLevel = Math.min(n - 1, 0);
                zoomMaxLevel = nZoomLevels - 1;
                //				zoomMaxLevel = Math.min(n+5, nZoomLevels-1);
                zoomCurrentLevel = n;
                break;
            }
        }

        zoomLevel = zoomLevels[zoomCurrentLevel];
        windowMdepthMax = StsMath.intervalRoundUp(mdepthMax, zoomLevel.unitsPerPixel);
        cursorY = displayHeight / 2;
        cursorMdepth = windowMdepthMax - cursorY*zoomLevel.unitsPerPixel;
    }

    public void relayoutWellWindow()
    {
        mdepthMin = well.getMinMDepth();
        mdepthMax = well.getMaxMDepth();
    }

    public int getCursorY()
    {
        return cursorY;
    }

    public double getMdepthFromMouseY(int mouseY, StsGLPanel glPanel)
    {
        return getMdepthFromGLMouseY(glPanel.getGLMouseY(mouseY));
    }

    public double getMdepthFromGLMouseY(int glMouseY)
    {
        return windowMdepthMax - glMouseY*zoomLevel.unitsPerPixel;
    }

    public double getCursorMdepth()
    {
        return cursorMdepth;
    }

    public String getZStringFromGLMouseY(int mouseY)
    {
        double mdepth = getMdepthFromGLMouseY(mouseY);
        return convertMDepthToZString(mdepth);
    }

    public String convertMDepthToZString(double mdepth)
    {
        String typeString = getUnits();
        String zLabel;
        double convertedZ;
        StsWellViewModel wellViewModel = this;

        if (typeString == wellViewModel.MDEPTH_TIME)
        {

            convertedZ = wellViewModel.well.getValueFromMDepth((float) mdepth, wellViewModel.MDEPTH);
            zLabel = wellViewModel.labelFormatter.format(convertedZ);

            convertedZ = wellViewModel.well.getValueFromMDepth((float) mdepth, wellViewModel.TIME);
            zLabel = zLabel + " " + wellViewModel.labelFormatter.format(convertedZ);

        }
        else if (typeString == wellViewModel.DEPTH_TIME)
        {

            convertedZ = wellViewModel.well.getValueFromMDepth((float) mdepth, wellViewModel.DEPTH);
            zLabel = wellViewModel.labelFormatter.format(convertedZ);

            convertedZ = wellViewModel.well.getValueFromMDepth((float) mdepth, wellViewModel.TIME);
            zLabel = zLabel + " " + wellViewModel.labelFormatter.format(convertedZ);

        }
        else
        {

            convertedZ = wellViewModel.well.getValueFromMDepth((float) mdepth, typeString);
            zLabel = wellViewModel.labelFormatter.format(convertedZ);

        }

        return zLabel;
    }

    public String getMdepthStringFromGLCursorY()
    {
        return getZStringFromGLMouseY(cursorY);
    }

    public void restoreZoom(int level)
    {
        zoomCurrentLevel = Math.min(level, zoomMaxLevel);
        zoomLevel = zoomLevels[zoomCurrentLevel];
    }

    public void zoomIn(StsView view, StsMouse mouse)
    {
        int glMouseY = view.getGLMouseY(mouse);
        mouseZ = windowMdepthMax - glMouseY * zoomLevel.unitsPerPixel;
        // System.out.println("zoomIn. currentLevel: " + zoomCurrentLevel);
        zoomCurrentLevel = Math.min(zoomCurrentLevel + 1, zoomMaxLevel);
        zoomLevel = zoomLevels[zoomCurrentLevel];
        windowMdepthMax = mouseZ + glMouseY * zoomLevel.unitsPerPixel;
        savePixels(false);
    }

    public void zoomOut(StsView view, StsMouse mouse)
    {
        int glMouseY = view.getGLMouseY(mouse);
        mouseZ = windowMdepthMax - glMouseY * zoomLevel.unitsPerPixel;
        // System.out.println("zoomOut. currentLevel: " + zoomCurrentLevel);
        zoomCurrentLevel = Math.max(zoomCurrentLevel - 1, zoomMinLevel);
        zoomLevel = zoomLevels[zoomCurrentLevel];
        windowMdepthMax = mouseZ + glMouseY * zoomLevel.unitsPerPixel;
        savePixels(false);
    }

    public void moveWindow(StsMouse mouse)
    {

        int dy = mouse.getMouseDelta().y;
        savePixels(false);
        if (dy == 0)
        {
            return;
        }
        windowMdepthMax += dy * zoomLevel.unitsPerPixel;
        savePixels(false);
        repaint();
    }

    public boolean setCursorPosition(int mouseY, StsWellView wellView)
    {
        moveCursor(mouseY, wellView);
        cursorPicked = true;
        return true;
    }

    public boolean checkCursorPicked(int mouseY, StsWellView wellView)
    {
        currentWellView = wellView;
        int newCursorY = wellView.getHeight() - mouseY; // distance up from bottom
        if (StsMath.near(newCursorY, cursorY, 2))
        {
            moveCursor(mouseY, wellView);
            cursorPicked = true;
        }
        else
        {
            cursorPicked = false;
        }
        return cursorPicked;
    }

    public void moveCursor(int mouseY)
    {
        moveCursor(mouseY, currentWellView);
    }

    public void moveCursor(int mouseY, StsWellView wellView)
    {
        currentWellView = wellView;
        cursorY = wellView.glPanel.getGLMouseY(mouseY);
        cursorMdepth = getMdepthFromGLMouseY(cursorY);
        repaint();
    }

    public boolean keyReleased(StsView view, StsMouse mouse, KeyEvent e)
    {
        return true;
    }

    public boolean keyPressed(StsView view, StsMouse mouse, KeyEvent e)
    {
        if (!mouse.isButtonDown(StsMouse.VIEW))
        {
            return false;
        }
        int keyCode = e.getKeyCode();
        if (keyCode == KeyEvent.VK_A || keyCode == KeyEvent.VK_Z || keyCode == KeyEvent.VK_UP)
        {
            zoomIn(view, mouse);
            cursorY = displayHeight / 2;
        }
        else if (keyCode == KeyEvent.VK_S || keyCode == KeyEvent.VK_X || keyCode == KeyEvent.VK_DOWN)
        {
            zoomOut(view, mouse);
        }
        repaint();
        return true;
    }

    public String[] getViewList()
    {
        // StsException.systemDebug(this, "getViewList", "called.");
        String[] operationList = new String[]{NO_ACTION};
        if (well.hasLogCurves())
        {
            operationList = (String[]) StsMath.arrayAddElement(operationList, ADD_LOG_CURVES);
            operationList = (String[]) StsMath.arrayAddElement(operationList, ADD_TRACE_DISPLAY);
            //operationList = (String[])StsMath.arrayAddElement(operationList, ADD_CURTAIN);
            operationList = (String[]) StsMath.arrayAddElement(operationList, ADD_SYNTHETIC);
            //operationList = (String[])StsMath.arrayAddElement(operationList, ADD_WELL_TRACK);
            operationList = (String[]) StsMath.arrayAddElement(operationList, ADD_TD_CURVE);
            operationList = (String[]) StsMath.arrayAddElement(operationList, ADD_INTERVAL_VELOCITY);
            operationList = (String[]) StsMath.arrayAddElement(operationList, EDIT_LOG_CURVE);
        }
        // if(currentModel.classHasObjects(StsSeismicVolume.class))
        // jbw     operationList = (String[])StsMath.arrayAddElement(operationList, ADD_CURTAIN );
        if (currentModel.stsClassHasThisFieldObject("com.Sts.DBTypes.StsVsp", "well", well))
            operationList = (String[]) StsMath.arrayAddElement(operationList, ADD_VSP);
        if (currentModel.stsClassHasThisFieldObject("com.Sts.DBTypes.StsTimeLogCurve", "well", well))
            operationList = (String[]) StsMath.arrayAddElement(operationList, ADD_DTS_PANEL);
        return operationList;
    }

    public StsVsp[] getAvailableVSPs()
    {
        if (!currentModel.classHasObjects(StsVsp.class)) return null;
        StsVsp[] vsps = (StsVsp[]) currentModel.getCastObjectList(StsVsp.class);
        int nVsps = vsps.length;
        StsVsp[] availVsps = new StsVsp[nVsps];
        int nn = 0;
        for (int n = 0; n < nVsps; n++)
            if (vsps[n].getWell() == well && !hasVsp(vsps[n])) availVsps[nn++] = vsps[n];
        if (nn == 0) return null;
        return (StsVsp[]) StsMath.trimArray(availVsps, nn);
    }

    private boolean hasVsp(StsVsp vsp)
    {
        for (int n = 0; n < wellViews.length; n++)
            return wellViews[n].displayedObject == vsp;
        return false;
    }

    public boolean hasVsps()
    {
        if (!currentModel.classHasObjects(StsVsp.class)) return false;
        StsVsp[] vsps = (StsVsp[]) currentModel.getCastObjectList(StsVsp.class);
        int nVsps = vsps.length;
        StsVsp[] availVsps = new StsVsp[nVsps];
        int nn = 0;
        for (int n = 0; n < nVsps; n++)
            if (vsps[n].getWell() == well && (!hasVsp(vsps[n]))) return true;
        return false;
    }

    public void setOperationType(String operationType)
    {

        if (operationType == ADD_LOG_CURVES)
        {
            selectLogCurves();
        }
        else if (operationType == ADD_SYNTHETIC)
        {
            String DT = selectLogCurvesForSynth("Sonic");
            if (DT == null) return;
            String RHOB = selectLogCurvesForSynth("Density");
            if (RHOB == null) return;
            StsLogCurvesView view = addLogTrackToWindow();
            view.setDisplaySynthetic(true);
            view.setSyntheticNames(DT, RHOB);
            view.addLogCurvePanel(0);
            view.track.logCurves.add(new StsLogCurve("SYNTHETIC"));
            // nSubWindows++;
            rebuild();
        }

        else if (operationType == ADD_TRACE_DISPLAY)
        {
            StsLogCurvesView view = addLogTrackToWindow();
            view.setDisplayTraces(true);
            view.addLogCurvePanel(0);
            view.track.logCurves.add(new StsLogCurve("SEISMIC"));
            // nSubWindows++;
            rebuild();
        }

        else if (operationType == ADD_TD_CURVE)
        {
            StsLogCurvesView view = addLogTrackToWindow();
            view.setDisplayDerived(true, "TIME/DEPTH");
            view.addLogCurvePanel(0);
            view.track.logCurves.add(new StsLogCurve("TIME/DEPTH"));
            // nSubWindows++;
            rebuild();
        }


        else if (operationType == ADD_INTERVAL_VELOCITY)
        {
            StsLogCurvesView view = addLogTrackToWindow();
            view.setDisplayDerived(true, "Interval Velocity");
            view.addLogCurvePanel(0);
            view.track.logCurves.add(new StsLogCurve("Interval Velocity"));
            // nSubWindows++;
            rebuild();
        }

        else if (operationType == EDIT_LOG_CURVE)
        {
            this.editLogCurves();
        }
        else if (operationType == ADD_WELL_TRACK)
        {
            addLogTrackToWindow();
            // nSubWindows++;
        }
        else if (operationType == ADD_VSP)
        {
            if (wellTrackView != null)
            {
                selectVSPs();
            }
        }
        else if (operationType == ADD_CURTAIN)
        {
            if (wellTrackView != null)
            {
                addSeismicCurtain();
                // nSubWindows++;
            }
        }
        else if (operationType == ADD_DTS_PANEL)
        {
            if (wellTrackView != null)
            {
                addDtsView();
                // nSubWindows++;
            }
            rebuild();
        }
    }

    public void display()
    {
        //seismicChanged();
        repaint();
    }

    public void display(StsGLPanel p)
    {
        display();
    }

    public void savePixels(boolean savem)
    {
        for (StsWellView wellView : wellViews)
            wellView.savePixels(savem);
    }

    public void repaintViews(Class viewClass)
    {
        if(wellViews == null) return;
        for(int n = 0; n < wellViews.length; n++)
        {
            if(viewClass.isAssignableFrom(wellViews[n].getClass()))
                wellViews[n].glPanel.repaint();
        }
    }

    public void repaint()
    {
        for (StsWellView wellView : wellViews)
            wellView.repaint();
    }

    public StsObject defaultEdit(StsMouse mouse)
    {
        return null;
    }

    public void setUnits(String units)
    {
        if ((units == TIME || (units == MDEPTH_TIME) || (units == DEPTH_TIME))
            && !currentModel.getProject().supportsTime())
        {
            StsMessageFiles.errorMessage("Project is not time capable, no td curves loaded for this well");
            return;
        }
        this.units = units;
    }

    public String getUnits()
    {
        return units;
    }

    static public Dimension getDefaultWellWindowSize()
    {
        return new Dimension(defaultPanelWidth, defaultHeight);
    }

    public void addVSPs(StsVsp[] vsps)
    {
        if (vsps == null) return;
        for (int i = 0; i < vsps.length; i++)
            addSelectedVSP(vsps[i]);
    }

    public void addSelectedVSP(StsVsp vsp)
    {
        StsVspView vspView = new StsVspView(this, currentModel, actionManager, nSubWindows, vsp);
        addWellView(vspView);
        rebuild();
    }

    //TODO redesign so that plug-in provides a viewConstructor with a static label "Add VSP Panel" which is included in drop-down list.
    //TODO When selected, this viewConstructor is called with this viewModel as an argument and performs actions below as needed on this wellViewModel.
    public void addSeismicCurtain()
    {
        StsSeismicCurtain seismicCurtain = well.getCreateSeismicCurtain();
        StsWellSeismicCurtainView seismicCurtainView = new StsWellSeismicCurtainView(this, currentModel, actionManager, nSubWindows, seismicCurtain);
        addWellView(seismicCurtainView);
    }

    //TODO as above
    public void addDtsView()
    {
        StsTimeLogCurve dataset = (StsTimeLogCurve) currentModel.getStsObjectWithThisFieldObject("com.Sts.DBTypes.StsTimeLogCurve", "well", well);
        if (dataset == null) return;
        StsTimeLogCurveWellView wellView = StsTimeLogCurveWellView.constructor(this, currentModel, actionManager, dataset, nSubWindows);
        addWellView(wellView);
    }

    public void removeWellWindowPanel(Component comp)
    {
        try
        {
            wellWindowPanel.remove(comp);
//            rebuild();
        }
        catch (Exception e)
        {
            StsException.outputException("StsLogCurvesView.removeLogCurve() failed.", e, StsException.WARNING);
        }
    }

    public void removeView(StsWellView wellView)
    {
        wellViews = (StsWellView[]) StsMath.arrayDeleteElement(wellViews, wellView);
		nSubWindows--;
    }

    /**
     * We might be interested in this object; returning true causes vewObjectChangedRepaint()
     * to be called which does the serious checking.
     */
    public boolean viewObjectChanged(Object source, Object object)
    {
        boolean changed = false;
        relayoutWellWindow();
        for (StsWellView wellView : wellViews)
            changed = changed | wellView.viewObjectChanged(source, object);

        // LiveWell - Adjust Zoom Position
        if(getWell() instanceof StsLiveWell)
        {
            if(currentModel.getProject().isRealtime())
                windowMdepthMax = StsMath.intervalRoundUp(mdepthMax, zoomLevel.unitsPerPixel);
        }
        return changed;
    }

    public boolean viewObjectRepaint(Object source, Object object)
    {
        boolean changed = false;
        for (StsWellView wellView : wellViews)
            changed = changed | wellView.viewObjectRepaint(source, object);
        return changed;
    }

    public void initAndStartGL()
    {
        for (StsWellView wellView : wellViews)
            wellView.initAndStartGL();
    }

    public StsModel getModel()
    {
        return currentModel;
    }

	public void pop2DOrtho(GL gl)
	{
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPopMatrix();
		gl.glMatrixMode(GL.GL_MODELVIEW);
	}

	public void push2DOrtho(GL gl, GLU glu, int width, int height)
	{
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPushMatrix();
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glLoadIdentity();
		glu.gluOrtho2D(0., width, 0., height);
		gl.glMatrixMode(GL.GL_MODELVIEW);
	}

    /** well measured depth max at bottom. */
    public double getWindowMdepthMax()
    {
        return windowMdepthMax;
    }

    /** well measured depth min at top of window. */
    public double getWindowMdepthMin(StsWellView wellView)
    {
        return windowMdepthMax - wellView.getHeight()*zoomLevel.unitsPerPixel;
    }
}
