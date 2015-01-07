
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.Collaboration.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.util.*;

public class StsCursor3dPanel extends StsJPanel implements ItemListener, ChangeListener //, StsSerializable
{
    StsModel model;
    private StsCursor3d cursor3d;
    StsRotatedGridBoundingBox rotatedBoundingBox;
    boolean isMainWindow = false;
    ButtonGroup sliderGroup = new ButtonGroup();
    boolean isEnabled = true;

    StsSliderBean sliderX = new StsSliderBean();
    StsSliderBean sliderY = new StsSliderBean();
    StsSliderBean sliderZ = new StsSliderBean();

    StsFloatFieldBean xBean = new StsFloatFieldBean(false, "X");
    StsFloatFieldBean yBean = new StsFloatFieldBean(false, "Y");

    int[] dirs = new int[] { XDIR, YDIR, ZDIR };
    public StsSliderBean[] sliders = new StsSliderBean[] { sliderX, sliderY, sliderZ };

    StsSliderBean currentSlider = null;
	int direction = NONE;
	boolean isDragging = false;
    private transient Vector changeListeners = null;
    private transient Vector itemListeners = null;
    public JCheckBox gridChk;
    public JCheckBox smallChk = new JCheckBox();
    StsWin3dFull window = null;

    static final boolean debug = false;

	static final int NONE = StsParameters.NO_MATCH;
	static final int XDIR = StsParameters.XDIR;
	static final int YDIR = StsParameters.YDIR;
	static final int ZDIR = StsParameters.ZDIR;

    /** object used to indicate that isGridCoordinateChanged checkbox has been switched.
     *  This object is passed by viewObjectChanged and viewObjectRepaint to views that have
     *  cursor axes displayed.
     */
    static public final String isGridCoordinateChanged = "gridCoordinateChanged";

    public StsCursor3dPanel()
    {
        super(false);
    }

    public StsCursor3dPanel(StsModel model, StsWin3dFull mainWindow, StsCursor3d cursor3d, Dimension size)
    {
        super(false);
        try
        {
            this.model = model;
			this.window = mainWindow;
            this.cursor3d = cursor3d;
            this.rotatedBoundingBox = model.getProject().getRotatedBoundingBox();
            this.isMainWindow = (mainWindow instanceof StsWin3d);
            constructPanel(size);
            setSliderValues();

            for (int n = 0; n < 3; n++)
                sliders[n].addToGroup(sliderGroup);
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

    protected void constructPanel(Dimension size) throws Exception
    {
//        setPreferredSize(size);
        sliderX.setIncrementLabel("Step");
        sliderX.setTextColor(Color.red);
        sliderX.setValueLabel("X");
 //       sliderX.setIcons("xSelected", "xDeselected");
        sliderX.setSelected(true);
        sliderY.setIncrementLabel("Step");
        sliderY.setTextColor(Color.green);
        sliderY.setValueLabel("Y");
        sliderY.setSelected(true);
 //       sliderY.setIcons("ySelected", "yDeselected");
        sliderZ.setIncrementLabel("Step");
        sliderZ.setTextColor(Color.blue);
        sliderZ.setValueLabel("Z");
        sliderZ.setSelected(true);

        smallChk.setText("Sliders only");
        smallChk.addItemListener(this);

        DecimalFormat format = new DecimalFormat("###0.0");
        xBean.setFormat(format);
        yBean.setFormat(format);
        xBean.setColumns(10);
        yBean.setColumns(10);
        xBean.gbc.anchor = GridBagConstraints.CENTER;
        yBean.gbc.anchor = GridBagConstraints.CENTER;

        StsJPanel checkPanel = StsJPanel.noInsets();
        checkPanel.gbc.anchor = GridBagConstraints.CENTER;
        checkPanel.addToRow(smallChk);
        if(isMainWindow)
        {
            gridChk = new JCheckBox();
            gridChk.setText("X-Y");
            gridChk.addItemListener(this);
            checkPanel.addEndRow(gridChk);
            boolean isCursorDisplayXY = window.getCursorDisplayXY();
            gridChk.setSelected(isCursorDisplayXY);
            int nSeismicVolumes = model.getNObjects(StsSeismicVolume.class);
			int nPreStackVolumes = model.getNObjects(StsPreStackLineSet3d.class);
            gridChk.setVisible(nSeismicVolumes + nPreStackVolumes > 0);
        }
        StsJPanel leftPanel =StsJPanel.noInsets();
        leftPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        leftPanel.add(checkPanel);
        leftPanel.add(xBean);
        leftPanel.add(yBean);

        gbc.fill = GridBagConstraints.NONE;
        addToRow(leftPanel);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        addToRow(sliderX);
        addToRow(sliderY);
        addEndRow(sliderZ);

        for (int n = 0; n < 3; n++)
        {
            sliders[n].addChangeListener(this);
            sliders[n].addItemListener(this);
        }
    }

    public synchronized void removeChangeListener(ChangeListener l)
    {
        if (changeListeners != null && changeListeners.contains(l))
        {
            Vector v = (Vector) changeListeners.clone();
            v.removeElement(l);
            changeListeners = v;
        }
    }

    public synchronized void addChangeListener(ChangeListener l)
    {
        Vector v = changeListeners == null ? new Vector(2) :
            (Vector) changeListeners.clone();
        if (!v.contains(l))
        {
            v.addElement(l);
            changeListeners = v;
        }
    }

    protected void fireStateChanged(ChangeEvent e)
    {
        if (changeListeners != null)
        {
            Vector listeners = changeListeners;
            int count = listeners.size();
            for (int i = 0; i < count; i++)
                ( (ChangeListener) listeners.elementAt(i)).stateChanged(e);
        }
    }

    public synchronized void removeItemListener(ItemListener l)
    {
        if (itemListeners != null && itemListeners.contains(l))
        {
            Vector v = (Vector) itemListeners.clone();
            v.removeElement(l);
            itemListeners = v;
        }
    }

    public synchronized void addItemListener(ItemListener l)
    {
        Vector v = itemListeners == null ? new Vector(2) :
            (Vector) itemListeners.clone();
        if (!v.contains(l))
        {
            v.addElement(l);
            itemListeners = v;
        }
    }

    protected void fireItemStateChanged(ItemEvent e)
    {
        if (itemListeners != null)
        {
            Vector listeners = itemListeners;
            int count = listeners.size();
            for (int i = 0; i < count; i++)
                ( (ItemListener) listeners.elementAt(i)).itemStateChanged(e);
        }
    }

	private void doSetGridCheckBoxVisibleAndSelected()
	{
        if(gridChk == null) return;
        if (gridChk.isVisible()) return;
        gridChk.setVisible(true);
        gridChk.setSelected(true);
        repaint();
	}

    public void setGridCheckboxState(boolean value_)
    {
        if(gridChk == null) return;
        final boolean value = value_;
        if(!model.getProject().rotatedBoundingBox.isRowColNumInitialized()) return;
        StsToolkit.runWaitOnEventThread(new Runnable() { public void run() { doSetGridCheckBoxSelected(value); }});
    }

    private void doSetGridCheckBoxSelected(boolean selected)
    {
        gridChk.setSelected(selected);
    }

    public void gridCheckBoxSetVisibleAndSelected()
    {
        if(!model.getProject().rotatedBoundingBox.isRowColNumInitialized()) return;
        StsToolkit.runWaitOnEventThread(new Runnable() { public void run() { doSetGridCheckBoxVisibleAndSelected(); }});
    }

    /** Programmatically set slider values from current cursor3d positions.
	 *  This assumes of course that cursor3d positions have already been set properly.
	 */

	public void setSliderValues()
    {
        StsToolkit.runLaterOnEventThread ( new Runnable() { public void run() { doSetValues(); } } );
//        if(rotatedBoundingBox.zMin == StsParameters.largeFloat || rotatedBoundingBox.zMax == -StsParameters.largeFloat) return;
//        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { doSetValues(); }});
    }

    public void doSetValues()
    {
        float cursorX = cursor3d.getCurrentDirCoordinate(StsCursor3d.XDIR);
        float cursorY = cursor3d.getCurrentDirCoordinate(StsCursor3d.YDIR);
        float cursorZ = cursor3d.getCurrentDirCoordinate(StsCursor3d.ZDIR);
        setSliderValues(cursorX, cursorY, cursorZ);
    }

    public void setSliderValues(float cursorX, float cursorY, float cursorZ)
    {
        boolean isGridCoordinates = !window.getCursorDisplayXY();

		// if we are changing coordinates, the change flips the currentDir, so reset it at bottom
		int currentDirNo = cursor3d.getCurrentDirNo();
        if(debug) StsException.systemDebug(this, "setSliderValues", " dirNo: " + currentDirNo +
                " cursorX: " + cursorX + " cursorY: " + cursorY + " cursorZ: " + cursorZ);
        if(!isGridCoordinates)
        {
            float[][] coorRanges = getCoordinateRanges();
            float value = getCoordinateValue(XDIR, cursorX);
            sliderX.initSliderValues(coorRanges[XDIR][0], coorRanges[XDIR][1], coorRanges[XDIR][2], value);
            value = getCoordinateValue(YDIR, cursorY);
            sliderY.initSliderValues(coorRanges[YDIR][0], coorRanges[YDIR][1], coorRanges[YDIR][2], value);
            value = getCoordinateValue(ZDIR, cursorZ);
            sliderZ.initSliderValues(coorRanges[ZDIR][0], coorRanges[ZDIR][1], coorRanges[ZDIR][2], value);
            sliderX.setValueLabel("X");
            sliderY.setValueLabel("Y");
        }
        else
        {
            float xLineMin = rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, rotatedBoundingBox.xMin);
            float yLineMin = rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, rotatedBoundingBox.yMin);
            float xLineMax = rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, rotatedBoundingBox.xMax);
            float yLineMax = rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, rotatedBoundingBox.yMax);

            sliderX.initSliderValues(xLineMin, xLineMax, rotatedBoundingBox.getColNumInc(),
                   rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, cursorX));
            sliderY.initSliderValues(yLineMin, yLineMax, rotatedBoundingBox.getRowNumInc(),
                   rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, cursorY));
            sliderZ.initSliderValues(rotatedBoundingBox.zMin, rotatedBoundingBox.zMax, rotatedBoundingBox.getZInc(),
                   rotatedBoundingBox.getNumFromCoor(StsCursor3d.ZDIR, cursorZ));

            sliderX.setValueLabel("XL");
            sliderY.setValueLabel("IL");
        }
		cursor3d.setCurrentDirNo(currentDirNo);
    }

    public float[][] getCoordinateRanges()
    {
        boolean isRotated = rotatedBoundingBox.getAngleSet() && rotatedBoundingBox.getAngle() != 0.0;
        float[][] minMaxs = new float[3][3];
        if(isRotated)
        {
            minMaxs[XDIR][0] =  rotatedBoundingBox.xMin;
            minMaxs[XDIR][1] =  rotatedBoundingBox.xMax;
            minMaxs[YDIR][0] =  rotatedBoundingBox.yMin;
            minMaxs[YDIR][1] =  rotatedBoundingBox.yMax;
        }
        else
        {
            double[] mins = rotatedBoundingBox.getAbsoluteXY(rotatedBoundingBox.xMin, rotatedBoundingBox.yMin);
            double[] maxs = rotatedBoundingBox.getAbsoluteXY(rotatedBoundingBox.xMax, rotatedBoundingBox.yMax);
            minMaxs[XDIR][0] =  (float)mins[0];
            minMaxs[XDIR][1] =  (float)maxs[0];
            minMaxs[YDIR][0] =  (float)mins[1];
            minMaxs[YDIR][1] =  (float)maxs[1];
        }
        minMaxs[XDIR][2] =  rotatedBoundingBox.xInc;
        minMaxs[YDIR][2] =  rotatedBoundingBox.yInc;
        minMaxs[ZDIR][0] =  rotatedBoundingBox.zMin;
        minMaxs[ZDIR][1] =  rotatedBoundingBox.zMax;
        minMaxs[ZDIR][2] =  rotatedBoundingBox.zInc;
        return minMaxs;
    }

    private float getCoordinateValue(int dir, float dirCoor)
    {
        float value;
        boolean isRotated = rotatedBoundingBox.getAngleSet() && rotatedBoundingBox.getAngle() != 0.0;
        if(isRotated) return dirCoor;

        if(dir == XDIR)
            return (float)model.getProject().getXOrigin() + dirCoor;
        else if(dir == YDIR)
            return (float)model.getProject().getYOrigin() + dirCoor;
        else
            return dirCoor;
    }

    public void step(int dir, boolean increment)
    {
        StsSliderBean slider = getSlider(dir);
        if(increment)
            slider.incrementSliderValue();
        else
            slider.decrementSliderValue();
    }

    private StsSliderBean getSlider(int dir)
    {
        if(dir == XDIR)
            return sliderX;
        else if(dir == YDIR)
            return sliderY;
        else
            return sliderZ;
    }

	public void setSliderValue(int dir, float dirCoor)
	{
        float[][] coorRanges;
        float value;
		boolean isGridCoordinates = !model.win3d.getCursorDisplayXY();

        coorRanges = getCoordinateRanges();
        value = getCoordinateValue(dir, dirCoor);
		if(!isGridCoordinates)
		{
			switch(dir)
			{
				case XDIR:
					sliderX.initSliderValues(coorRanges[XDIR][0], coorRanges[XDIR][1], coorRanges[XDIR][2], value);
					sliderX.setValueLabel("X");
					break;
				case YDIR:
					sliderY.initSliderValues(coorRanges[YDIR][0], coorRanges[YDIR][1], coorRanges[YDIR][2], value);
					sliderY.setValueLabel("Y");
					break;
				case ZDIR:
					sliderZ.initSliderValues(coorRanges[ZDIR][0], coorRanges[ZDIR][1], coorRanges[ZDIR][2], value);
			}
		}
		else
		{
			switch(dir)
			{
				case XDIR:
					float xLineMin = rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, rotatedBoundingBox.xMin);
					float xLineMax = rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, rotatedBoundingBox.xMax);
					sliderX.initSliderValues(xLineMin, xLineMax, rotatedBoundingBox.getColNumInc(), rotatedBoundingBox.getNumFromCoor(StsCursor3d.XDIR, dirCoor));
					sliderX.setValueLabel("XL");
					break;
				case YDIR:
					float yLineMin = rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, rotatedBoundingBox.yMin);
					float yLineMax = rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, rotatedBoundingBox.yMax);
					sliderY.initSliderValues(yLineMin, yLineMax, rotatedBoundingBox.getRowNumInc(),
						   rotatedBoundingBox.getNumFromCoor(StsCursor3d.YDIR, dirCoor));
                    sliderY.setValueLabel("IL");
					break;
				case ZDIR:
					sliderZ.initSliderValues(rotatedBoundingBox.zMin, rotatedBoundingBox.zMax, rotatedBoundingBox.getZInc(),
						   rotatedBoundingBox.getNumFromCoor(StsCursor3d.ZDIR, dirCoor));
			}
		}
    }

    public StsSliderBean getSliderX() { return sliderX; }
    public StsSliderBean getSliderY() { return sliderY; }
    public StsSliderBean getSliderZ() { return sliderZ; }

    public void selectSlider(StsSliderBean slider)
    {
        if( slider == currentSlider ) return;
        currentSlider = slider;
        int direction = NONE;
        for( int i=0; i<dirs.length; i++ )
        {
            if( slider == sliders[i] )
            {
                direction = dirs[i];
                break;
            }
        }
//        slider.setButtonSelected(true);
        ItemEvent event = new ItemEvent(slider.getCheckBoxSlider(), direction, slider, ItemEvent.SELECTED);
        fireItemStateChanged(event);
    }

    public void itemStateChanged(ItemEvent e)
    {
        ItemEvent event = null;
        int stateChange = e.getStateChange();
        boolean selected = (stateChange == ItemEvent.SELECTED);
        Object source = e.getSource();

        if(source == gridChk)
        {
//            boolean isSelected = gridChk.isSelected();
            boolean cursorDisplayXY = window.getCursorDisplayXY();
            if(cursorDisplayXY == selected) return;
            window.setCursorDisplayXY(selected);
            model.resetAllSliderValues();
            model.viewObjectRepaint(this, isGridCoordinateChanged);
        }
        else if(source == smallChk)
        {
            sliderX.sliderOnly(selected);
            sliderY.sliderOnly(selected);
            sliderZ.sliderOnly(selected);
            if(selected)
            {
                gridChk.setVisible(false);
                xBean.setVisible(false);
                yBean.setVisible(false);
            }
            else
            {
                gridChk.setVisible(true);
                xBean.setVisible(true);
                yBean.setVisible(true);
            }
            window.validate();
            window.repaint();
        }
        else if(source instanceof StsSliderBean)
        {
            StsSliderBean slider = (StsSliderBean)source;
            Object item = e.getItem();
            if (item instanceof JCheckBox)
                checkBoxSelected(slider, selected);
            else if (item instanceof JToggleButton)
                buttonSelected(slider, selected);
        }
    }

    public void buttonSelected(StsSliderBean slider, boolean enable)
    {
        slider.setButtonSelected(enable);
        if(!enable) return;
        setSliderSelected(slider, enable);
        if(!slider.isSelected()) slider.setCheckBoxModelSelected(true);
        win3dDisplayAll();
    }

    public void checkBoxSelected(StsSliderBean slider, boolean enable)
    {
        if (enable)
        {
            if (currentSlider != slider)
            {
                setSliderSelected(slider, true);
                win3dDisplayAll();
//                event = new ItemEvent(e.getItemSelectable(), direction, slider, stateChange);
//                fireItemStateChanged(event);
            }
        }
        else
        {
            setSliderSelected(slider, false);
            if (currentSlider == slider)
            {
                currentSlider = null;
                for (int n = 0; n < 3; n++)
                {
                    if (sliders[n] != slider && sliders[n].isSelected())
                    {
                        setSliderSelected(sliders[n], true);
                        break;
                    }
                }
            }
            win3dDisplayAll();
        }
    }

    private void win3dDisplayAll()
    {
        StsWindowFamily family = model.getWindowFamily(window);
        family.win3dDisplay();
    }

	public void setSliderSelected(int i, boolean enable)
   {
       StsWindowFamily family = model.getWindowFamily(window);
       family.setSelectedDirection(i, enable);
       sliders[i].setCheckBoxSelected(enable);
    }

    public void setSliderSelected(StsSliderBean slider, boolean enable)
    {
        if(enable)
        {
            if (currentSlider == slider) return;
            currentSlider = slider;
        }
		setSliderState(slider);
        StsWindowFamily family = model.getWindowFamily(window);
        family.setSelectedDirection(direction, enable);
        if(debug) System.out.println("Direction " + StsParameters.coorLabels[direction] + " enabled: " + enable);
    }

	private void setSliderState(StsSliderBean slider)
	{
		direction = getSliderDirection(slider);
		// We could be dragging a cursor3d slider or selecting a point in the 3d view
		isDragging = slider.isAdjusting || window.hasFocus() && StsGLPanel.isLeftMouseDown();
		if(debug) StsException.systemDebug(this, "setSliderState", "stateChanged for slider " + direction + " slider dragging: " +
                slider.isAdjusting + " window has focus " + window.hasFocus() + " leftMouseDown " + StsGLPanel.isLeftMouseDown());
		cursor3d.setSliderState(direction, isDragging);

        double[] xy = {sliders[0].getValue(), sliders[1].getValue()};
        if(rotatedBoundingBox.getAngleSet())
        	xy = model.getProject().getAbsoluteXYCoordinates(new float[] { sliders[0].getValue(), sliders[1].getValue() });
        xBean.setValue(xy[0]);
        yBean.setValue(xy[1]);
	}

    private int getSliderDirection(StsSliderBean slider)
    {
        for(int n = 0; n < 3; n++)
            if(sliders[n] == slider) return n;
        return StsParameters.NONE;
    }

    public StsWin3dBase getWindow()
    {
        return window;
    }
/*
    public void stateChanged(ChangeEvent e)
    {
        if((e.getSource() instanceof StsSliderBean) )
        {
            StsSliderBean slider = (StsSliderBean) e.getSource();
            setSliderSelected(slider, true);
//            selectSlider(slider);
            fireStateChanged(e);
        }
        else
            ;
    }
 */
    public void stateChanged(ChangeEvent e)
    {
        if( !(e.getSource() instanceof StsSliderBean) ) return;

        StsSliderBean slider = (StsSliderBean) e.getSource();
		setSliderState(slider);
        if( slider.adjustingChanged() )
        {
            StsWindowFamily windowFamily = window.getWindowFamily();
            Iterator viewIterator = windowFamily.getWindowViewIteratorOfType(StsView3d.class);
            while(viewIterator.hasNext())
            {
                StsView3d view3d = (StsView3d)viewIterator.next();
                view3d.set3dOverlay(isDragging);
            }
        }
//        if(isDragging)
        {
			float sliderValue = getSliderValue(slider, direction);
            setCursor(direction, sliderValue);
			if(debug) StsException.systemDebug(this, "setCursor", " direction: " + direction + " value: " + sliderValue);
            window.getWindowFamily().win3dDisplay();
        }
    }

    public void lockSlider(int direction)
    {
        sliders[direction].lock();
    }

    public void unlockSlider(int direction)
    {
        sliders[direction].unlock();
    }

    private void setCursor(int dir, float dirCoor)
	{
		StsWindowFamily family = model.getWindowFamily(window);
		if (family == null) return;
        cursor3d.adjustCursor(dir, dirCoor);
        collaborationSetCursor(family, dir, dirCoor);
    }

    private void collaborationSetCursor(StsWindowFamily family, int dir, float dirCoor)
    {
        StsWin3dBase parent = family.getParent();
        StsCollaboration collaboration = StsCollaboration.getCollaboration();
		if (parent == window && collaboration != null && collaboration.hasPeers())
			model.createTransientCursorChangeTransaction(dir, dirCoor);
	}

    private float getSliderValue(StsSliderBean slider, int currentDirection)
    {
        float value = slider.getValue();
        boolean isRotated = false;
        if((rotatedBoundingBox.getAngleSet()) && (rotatedBoundingBox.getAngle() != 0.0))
            isRotated = true;
        if(cursor3d.getIsGridCoordinates())
        {
            if (currentDirection == StsCursor3d.XDIR)
                value = rotatedBoundingBox.getXFromColNum(value);
            else if (currentDirection == StsCursor3d.YDIR)
                value = rotatedBoundingBox.getYFromRowNum(value);
        }
        else
        {
            if(!isRotated)
            {
                // convert to local coordinates
                if(currentDirection == XDIR)
                    value = value - (float)model.getProject().getXOrigin();
                else if(currentDirection == YDIR)
                    value = value - (float)model.getProject().getYOrigin();
            }
        }
        return value;
    }
    /*
    private float getSliderValue(StsSliderBean slider, int currentDirection)
    {
//				System.out.println("Value changed: " + (float)slider.getValue());
//                if(currentDirection == -1) return;
        float sliderValue = computeRelativeSliderValue(slider, currentDirection);

        if(cursor3d.getIsGridCoordinates() && (currentDirection == StsCursor3d.XDIR || currentDirection == StsCursor3d.YDIR) )
        {
			if(rotatedBoundingBox != null)
			{
                 if (currentDirection == StsCursor3d.XDIR)
                    sliderValue = rotatedBoundingBox.getXFromColNum(sliderValue);
                else if (currentDirection == StsCursor3d.YDIR)
                    sliderValue = rotatedBoundingBox.getYFromRowNum(sliderValue);
            }
        }
        return sliderValue;
    }
     */
    public void setCurrentDirNo(int dirNo)
    {
        switch(dirNo)
        {
            case StsCursor3d.XDIR:
                sliderX.setButtonModelSelected(true);
                break;
            case StsCursor3d.YDIR:
                sliderY.setButtonModelSelected(true);
                break;
            case StsCursor3d.ZDIR:
                sliderZ.setButtonModelSelected(true);
                break;
        }
    }

	public void setEditable(boolean enabled)
	{
        if(isEnabled == enabled) return;
        isEnabled = enabled;

        Runnable runnable = new Runnable()
        {
            public void run()
            {
                sliderX.setEditable(isEnabled);
                sliderY.setEditable(isEnabled);
                sliderZ.setEditable(isEnabled);
            }
		};
        StsToolkit.runLaterOnEventThread(runnable);
    }

	static final long serialVersionUID = 1l;

    public StsCursor3d getCursor3d()
    {
        return cursor3d;
    }

    public void setCursor3d(StsCursor3d cursor3d)
    {
        this.cursor3d = cursor3d;
    }
}