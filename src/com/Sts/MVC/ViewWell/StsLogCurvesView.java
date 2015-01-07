package com.Sts.MVC.ViewWell;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.Actions.StsAction;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;
import java.util.ArrayList;

public class StsLogCurvesView extends StsWellView
{
    public StsLogTrack track;
    private boolean displayTraces = false;
    private boolean displaySynthetic = false;
    private boolean displayDerived = false;
    private String derivedLog = null;
    transient boolean displayValues = true;
    private String DTname = null;
    private String RHOBname = null;
    private transient JPopupMenu popup = null;
    ArrayList<StsPoint> drawPoints = null;

    transient public StsLogCurveNamePanel[] curveNamePanels;

    static public final String viewLogTrack = "Well Log Track View";

    public StsLogCurvesView()
    {

    }

    public StsLogCurvesView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow)
    {
        this(wellViewModel, model, actionManager, nSubWindow, new StsLogTrack(wellViewModel, model, nSubWindow));
    }

    public StsLogCurvesView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow, StsLogTrack track)
    {
        super(track);

        this.track = track;
        initializeView(wellViewModel, model, actionManager, nSubWindow);
    }

    public boolean initializeView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow)
    {

        try
        {
            this.wellViewModel = wellViewModel;
            this.model = model;
            this.well = wellViewModel.well;

            if(this.track != null)
            {
                track.setDisplayTraces(this.displayTraces);
                track.setDisplaySynthetic(this.displaySynthetic);
                track.setDisplayDerived(this.displayDerived, this.derivedLog);
                track.setSyntheticNames(DTname, RHOBname);
            }
            StsWellWindowPanel wellWindowPanel = wellViewModel.getWellWindowPanel();
            curveNameBackPanel = new StsJPanel();
            curveNameBackPanel.setBorder(BorderFactory.createRaisedBevelBorder());

            innerPanel = wellWindowPanel.getNewPanel(wellViewModel.curveTrackWidth, wellViewModel.displayHeight);

			GridBagConstraints gbc = innerPanel.gbc;
            gbc.weightx = 0.1;
            gbc.weighty = 0;
            gbc.gridy = 0;
            gbc.gridx = nSubWindow;
            gbc.fill = GridBagConstraints.HORIZONTAL;
            curveNameBackPanel.setMaximumSize(new Dimension(1000, 40));
            // wellWindowPanel.innerPanel.add(curveNameBackPanel);
            innerPanel.add(curveNameBackPanel);
            curveNamePanels = new StsLogCurveNamePanel[StsLogTrack.maxCurvesPerTrack];
            if(Main.useJPanel)
                glPanel = new StsGLJPanel(model, actionManager, wellViewModel.curveTrackWidth, wellViewModel.displayHeight, this);
            else
                glPanel = new StsGLPanel(model, actionManager, wellViewModel.curveTrackWidth, wellViewModel.displayHeight, this);
            // jbw glPanel.setBorder(BorderFactory.createLoweredBevelBorder());
            // glPanel.setBackgroundColor(Color.WHITE);
            gbc.gridy = 2;
            gbc.weightx = 1;
            gbc.weighty = 1;
            gbc.fill = GridBagConstraints.BOTH;
            gbc.anchor = GridBagConstraints.LAST_LINE_START;
            innerPanel.add(glPanel);
        //GraphicsDevice gDevice = glPanel.getGraphicsConfiguration().getDevice();
        	glPanel.initAndStartGL(glPanel.getGraphicsDevice());
            gbc.anchor = GridBagConstraints.CENTER;
            rebuildView(wellViewModel, model, nSubWindow);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLogCurvesView.constructor() failed.", e, StsException.WARNING);
            return false;
        }
    }

    public void toggleSeismic()
    {
        track.displaySeismic = !track.displaySeismic;
    }

    public void setDisplayTraces(boolean b)
    {
        this.displayTraces = b;
        track.displayTraces = b;
    }

    public void setDisplaySynthetic(boolean b)
    {
        this.displaySynthetic = b;
        track.displaySynthetic = b;
    }

    public void setSyntheticNames(String DT, String RHOB)
    {
        this.DTname = DT;
        this.RHOBname = RHOB;
        track.setSyntheticNames(DT, RHOB);
    }

    public void setDisplayDerived(boolean b, String Derived)
    {
        this.displayDerived = b;
        this.derivedLog = Derived;
        track.displayDerived = b;
        track.derivedLog = Derived;
    }

    public StsLogTrack getTrack()
    {
        return track;
    }

    public void reshape(int x, int y, int width, int height)
    {
        wellViewModel.setCurveTrackSize(new Dimension(width, height));
        viewPortChanged(x, y, width, height);
    }

    public void display(GLAutoDrawable component)
    {
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }

        //TODO Figure out how we can get the size change event regisered */
        if(glPanel.viewChanged)
        {
            computeProjectionMatrix();
            glPanel.viewChanged = false;
        }
        clearToBackground(GL.GL_COLOR_BUFFER_BIT);
        updateCurveNamePanels();
        track.drawLogTrack(component, displayValues, gl, glu, glPanel.getMousePoint(), getWidth(), getHeight(), wellViewModel, this, drawPoints);
//        drawMarkers(gl,glu, 10);
        drawCursor(gl, 0, glPanel.getGLWidth(), glPanel.getGLHeight(), false);
    }

    private void updateCurveNamePanels()
    {
        for(int n = 0; n < curveNamePanels.length; n++)
            if(curveNamePanels[n] != null) curveNamePanels[n].setValueLabels();
    }

    public boolean viewObjectChanged(Object source, Object object)
	{
        glPanel.repaint();
		return true;
	}
	/** object being viewed is changed. Repaint this view if affected.
	 * Implement as needed in concrete subclasses.
	 */
	public boolean viewObjectRepaint(Object source, Object object)
	{
        glPanel.repaint();
		return true;
	}

    public void setAxisRanges()
    {

    }

    public void setDefaultView()
    {
    }

    public void addLogs()
    {
        wellViewModel.selectLogCurves(this);
    }

    public void removeLogs()
    {
        wellViewModel.removeLogPanel(this);
    }

    public boolean addLogCurve(StsLogCurve logCurve)
    {
        if(logCurve == null) return false;
        int nCurves = track.logCurves.getSize();
        track.logCurves.add(logCurve);
        addLogCurve(logCurve, nCurves);
        track.initializeSelectedCurve();
        return true;
    }

    public int getNLogCurves()
    {
        return track.logCurves.getSize();
    }

    public void addLogCurve(StsLogCurve logCurve, int nCurve)
    {
        try
        {
            StsLogCurveNamePanel curveNamePanel = new StsLogCurveNamePanel(logCurve);
            if(curveNamePanel == null) return;
            track.addRadioButtonToGroup(curveNamePanel.selectCurveRadioButton);
            curveNamePanels[nCurve] = curveNamePanel;
            curveNameBackPanel.gbc.fill = GridBagConstraints.BOTH;
            curveNameBackPanel.add(curveNamePanel);
        }
        catch(Exception e)
        {
            StsException.outputException("StsLogCurvesView.addLogCurve() failed.", e, StsException.WARNING);
        }
    }

    public void addLogCurvePanel(int nCurve)
    {
        try
        {
            StsLogCurveNamePanel curveNamePanel = constructCurveNamePanel();
            if(curveNamePanel == null) return;
            curveNamePanels[nCurve] = curveNamePanel;
            curveNameBackPanel.gbc.fill = GridBagConstraints.BOTH;
            curveNameBackPanel.add(curveNamePanel);
        }
        catch(Exception e)
        {
            StsException.outputException("StsLogCurvesView.addLogCurve() failed.", e, StsException.WARNING);
        }
    }

    public boolean rebuildView(StsWellViewModel wellViewModel, StsModel model, int nSubWindow)
    {
        try
        {
            this.wellViewModel = wellViewModel;
            this.model = model;

            addLogCurvesToView();
            return true;
        }
        catch(Exception e)
        {
            StsException.systemError("StsLogTrack.rebuildView() failed.");
            return false;
        }
    }

    private void addLogCurvesToView()
    {
        int nCurves = track.logCurves.getSize();
        for(int n = 0; n < nCurves; n++)
        {
            StsLogCurve logCurve = (StsLogCurve) track.logCurves.getElement(n);
            addLogCurve(logCurve, n);
        }
        track.initializeSelectedCurve();
    }

    public void removeLogCurve(StsLogCurve logCurve)
    {
        try
        {
            int removeIdx = 0;

            // Remove log from panel.
            for(int i = 0; i < curveNamePanels.length; i++)
            {
                if(curveNamePanels[i] != null)
                {
                    if(curveNamePanels[i].getName().equals(logCurve.getName()))
                    {
                        removeIdx = -1;
                        curveNameBackPanel.remove(curveNamePanels[i]);
                    }
                    else
                        curveNamePanels[i + removeIdx] = curveNamePanels[i];
                }
            }
            if(removeIdx == -1)
                curveNamePanels[curveNamePanels.length - 1] = null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLogCurvesView.removeLogCurve() failed.", e, StsException.WARNING);
        }
        if(track.removeLogCurve(logCurve))
            removePanel();
    }

    private StsLogCurveNamePanel constructCurveNamePanel()
    {
        String curveName = "";
        float min = 0;
        float max = 1;
        StsPoint[] points = wellViewModel.well.getRotatedPoints();
        int nPoints = points.length;
        float t0 = points[0].getT();
        float m0 = points[0].getM();
        float t1 = points[nPoints - 1].getT();
        float m1 = points[nPoints - 1].getM();

        t1 = t1 * 1.1f;
        if(t0 < 0) t0 = 0;
        // assume monotonic;

        if(displayTraces)
            curveName = "Seismic";
        if(displaySynthetic)
            curveName = "Synthetic Seismic";
        if(displayDerived)
        {
            curveName = derivedLog;
            if(derivedLog.equals("TIME/DEPTH"))
            {
                min = t0;
                max = t1;
            }
            else
            {
                min = 000.f;
                max = 20000.f;
            }
        }
        StsLogCurveNamePanel curveNamePanel = new StsLogCurveNamePanel(curveName, min, max, "LIN");
        track.addRadioButtonToGroup(curveNamePanel.selectCurveRadioButton);
        return curveNamePanel;
    }

    private StsJPanel constructCurveNamePanel(StsLogCurve logCurve)
    {
        if(displayDerived || displayTraces || displaySynthetic) return constructCurveNamePanel();

        StsLogCurveType logCurveType = logCurve.getLogCurveType();
        String curveName = logCurve.getName();
        StsJPanel curveNamePanel = new StsJPanel();
        curveNamePanel.setName(curveName);
        curveNamePanel.setBorder(BorderFactory.createRaisedBevelBorder());
        //        curveNamePanel.setMinimumSize(new Dimension(curveTrackWidth, 25));
        //        curveNamePanel.setPreferredSize(new Dimension(curveTrackWidth, 25));
        JRadioButton selectCurveRadioButton = new JRadioButton();
        selectCurveRadioButton.setMargin(new Insets(0, 0, 0, 2));
//        selectCurveRadioButton.setPreferredSize(new Dimension(15, 15));
//        selectCurveRadioButton.setMinimumSize(new Dimension(15, 15));
        selectCurveRadioButton.setContentAreaFilled(false);
        track.addRadioButtonToGroup(selectCurveRadioButton);
        selectCurveRadioButton.setActionCommand(curveName);
        JLabel minValueLabel = new JLabel();
        JLabel curveLabel = new JLabel();
        JLabel maxValueLabel = new JLabel();
        JLabel scaleTypeLabel = new JLabel();
        curveLabel.setText(curveName);
        if(logCurveType != null)
        {
            String minValueString = StsMath.formatNumber(logCurveType.getScaleMin(), 5, 5);
            minValueLabel.setText(minValueString);
            String maxValueString = StsMath.formatNumber(logCurveType.getScaleMax(), 5, 5);
            maxValueLabel.setText(maxValueString);
            scaleTypeLabel.setText(logCurveType.getScaleTypeString());
            StsColor color = logCurve.getStsColor();
            curveNamePanel.setBackground(color.getColor());
            selectCurveRadioButton.setForeground(color.getColor());
        }

        curveNamePanel.addToRow(selectCurveRadioButton);
        curveNamePanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        curveNamePanel.addToRow(minValueLabel);
        curveNamePanel.gbc.fill = GridBagConstraints.NONE;
        curveNamePanel.addToRow(curveLabel);
        curveNamePanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        curveNamePanel.addToRow(maxValueLabel);
        curveNamePanel.gbc.fill = GridBagConstraints.NONE;
        curveNamePanel.addToRow(scaleTypeLabel);
        return curveNamePanel;
    }

    public void removePanel()
    {
        curveNameBackPanel.removeAll();
        StsWellWindowPanel wellWindowPanel = wellViewModel.getWellWindowPanel();
        //wellWindowPanel.removeInner(innerPanel);
        innerPanel.remove(curveNameBackPanel);
        glPanel.removeAll();

        wellViewModel.removeWellWindowPanel(glPanel);
        wellViewModel.removeView(this);
        wellWindowPanel.removeInner(innerPanel);
		wellWindowPanel.rebuild();
        //innerPanel.setVisible(false);
        this.delete();
        wellViewModel.rebuild();
    }

    public void viewChangedRepaint()
    {
        if(glPanel != null)
            glPanel.repaint();
        else
            System.out.println("null");
    }

    public StsLogTrack getLogTrack()
    {
        return track;
    }

    public void delete()
    {
        curveNameBackPanel = null;
        curveNamePanels = null;
        glPanel = null;
        track.delete();
    }

    public void savePixels(boolean b)
    {
        glPanel.savePixels(b);
    }

    public String getViewClassname()
    {
        return viewLogTrack;
    }

    public void resetToOrigin()
    {
    }

    public byte getHorizontalAxisType()
    {
        return AXIS_TYPE_NONE;
    }

    public byte getVerticalAxisType()
    {
        return AXIS_TYPE_MEASURED_DEPTH;
    }

    public void showPopupMenu(StsMouse mouse)
    {
        popup = new JPopupMenu();
        glPanel.add(popup);

        StsMenuItem toggleSeismicDisplay = new StsMenuItem();
        StsMenuItem removeLogs = new StsMenuItem();
        StsMenuItem addLogs = new StsMenuItem();
        StsMenuItem removePanel = new StsMenuItem();

        toggleSeismicDisplay.setMenuActionListener("Toggle Seismic", this, "toggleSeismic", null);
        removeLogs.setMenuActionListener("Remove Logs", this, "removeLogs", null);
        addLogs.setMenuActionListener("Add Logs", this, "addLogs", null);
        removePanel.setMenuActionListener("Remove Panel", this, "removePanel", null);

        popup.add(toggleSeismicDisplay);
        popup.add(removeLogs);
        popup.add(addLogs);
        popup.add(removePanel);

        popup.show(glPanel, mouse.getX(), mouse.getY());
    }

    public void cancelPopupMenu(StsMouse mouse)
    {
        if(popup != null) popup.setVisible(false);
    }

    public void performMouseAction(StsActionManager actionManager, StsMouse mouse)
    {
        try
        {
            StsAction currentAction = actionManager.getCurrentAction();
            if(currentAction != null && currentAction.performMouseAction(mouse, this))
                return;

            int currentButton = mouse.getCurrentButton();
            if(currentButton == StsMouse.LEFT)
            {
                int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
                if(buttonState == StsMouse.PRESSED)
                {
                    displayValues = true;
                    if(!wellViewModel.cursorPicked)
                        wellViewModel.checkCursorPicked(wellViewModel.getCursorY(), this); // jbw
                }
                else if(buttonState == StsMouse.DRAGGED)
                {
                    displayValues = true;
                    if(wellViewModel.cursorPicked)
                        wellViewModel.moveCursor(mouse.getMousePoint().y, this);
                    else
                    {
                        // wellViewModel.setCursorPanel(glPanel);
                        glPanel.setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
                        viewChangedRepaint();
                    }
                }
                else if(buttonState == StsMouse.RELEASED)
                {
                    glPanel.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                    wellViewModel.cursorPicked = false;
                    displayValues = false;
                    savePixels(false);
                    wellViewModel.display();
                }
            }

            // If any right mouse action, move view
            else if(currentButton == StsMouse.RIGHT)
            {
                int buttonState = mouse.getButtonStateCheckClear(StsMouse.RIGHT);
                if((buttonState == StsMouse.DRAGGED))
                    wellViewModel.moveWindow(mouse);
            }

            // If middle mouse button clicked, terminate any active function.
            else if (currentButton == StsMouse.MIDDLE)
            {
                if(mouse.getButtonStateCheckClear(StsMouse.POPUP) == StsMouse.PRESSED)
                {
                    showPopupMenu(mouse);
                    //mouse.clearButtonState(StsMouse.MIDDLE, StsMouse.CLEARED);
                }
                if(mouse.getButtonStateCheckClear(StsMouse.POPUP) == StsMouse.RELEASED)
                {
                    cancelPopupMenu(mouse);
                    //mouse.clearButtonState(StsMouse.MIDDLE, StsMouse.CLEARED);
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "performMouseAction", e);
        }
    }

    public void setDrawPoints(ArrayList<StsPoint> drawPoints)
    {
        this.drawPoints = drawPoints;
    }
}
