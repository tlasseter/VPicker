
//Title:        S2S: Seismic-to-simulation Slider
//Version:
//Copyright:    Copyright (c) 2001
//Author:       T J Lasseter
//Company:      4D Systems LLC
//Description:  Slider with value and increment controls


package com.Sts.UI;

import com.Sts.Utilities.StsException;
import com.Sts.Utilities.StsMath;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.io.*;

public class StsJSlider extends JSlider
{
    private int increment = 0;
    private int lastValue = 0;
    private boolean isDragging = false;
    /** Slider can be programmatically locked; uses LockableModelListener inner class */
    private boolean locked = false;
    boolean debug = false;
    String debugName;

    public StsJSlider()
    {
/*
		try
		{
			SliderUI sliderUI = this.getUI();
			if(sliderUI instanceof BasicSliderUI)
			{
				BasicSliderUI basicSliderUI = (BasicSliderUI)sliderUI;
				basicSliderUI.scrollListener.setScrollByBlock(false);
			}
			Class c = sliderUI.getClass();
			Method method = sliderUI.getClass().getMethod("setScrollByBlock", new Class[] { Boolean.TYPE });
	        if(method == null) return;
			method.invoke(sliderUI, new Object[] { new Boolean(false) });
		}
		catch(Exception e)
		{
			StsException.outputException("StsJSlider.constructor() failed.", e, StsException.WARNING);
		}
*/
	}

    public void setDebug()
    {
        debug = true;
    }

    public void setDebug(String sliderName)
    {
        this.debugName = sliderName;
    }

    public void initSliderValues(int min, int max, int increment, int value)
    {
        if(debug) StsException.systemDebug(this, "initSliderValues", " slider " + debugName + " value changed from " + lastValue + " to " + value);
        snapToTicks = true;
        BoundedRangeModel m = getModel();
        m.setRangeProperties(value, 0, min, max, false);
		boolean drawTicks = increment != 0 && (max - min)/increment < 200;
/*
		while(((max-min)/increment) > 5000)
		{
			increment *= 2;
		}
*/
		setIncrement(increment, drawTicks);
        lastValue = value;
    }
/*
    public void setValue(int value)
    {
    	isDragging = getValueIsAdjusting();
		if(value == lastValue) return;
		int change = Math.abs(value - lastValue);
        lastValue = value;
        super.setValue(value);
    }
*/
	public void setValue(int value)
	 {
		 if(value == lastValue) return;
         if(debug) StsException.systemDebug(this, "setValue", " slider " + debugName + " value changed from " + lastValue + " to " + value);
		 isDragging = getValueIsAdjusting();
		 if( increment != 0 )
		 {
			 // try to compensate for the bug in the block scroll and snapToTicks
			 // no matter what the tick settings are it always scrolls 10%
             int max = getMaximum();
             int min = getMinimum();
			 int blockIncrement = (max - min)/10;
			 int change = value - lastValue;
			 if(blockIncrement == change) // round up
                 value = StsMath.minMax(StsMath.intervalRoundUpAbove(lastValue, increment), min, max);
             else if(blockIncrement == -change)
                 value = StsMath.minMax(StsMath.intervalRoundDownBelow(lastValue, increment), min, max);
			 lastValue = value;
		 }
		 super.setValue(value);
	 }

    public void computeNewValue()
    {
		setValue(super.getValue());
    }

    public void setIncrement(int newIncrement, boolean drawTicks)
    {
//		System.out.println("StsJSlider.setIncrement called with increment " + newIncrement + " on thread " + Thread.currentThread().getName());
        increment = newIncrement;
		if(drawTicks)
            setMajorTickSpacing(increment);
        else
            setMajorTickSpacing(0); // don't draw major ticks if too many
        setMinorTickSpacing(0); // never draw minor ticks
	}

    public int getIncrement()
    {
        return increment;
    }

    public boolean isDraggingSlider()
    {
        return getValueIsAdjusting();
//        return isDragging;
    }

    public void lock() { locked = true; }
    public void unlock() { locked = false; }

    protected ChangeListener createChangeListener()
    {
        return new LockableModelListener();
    }

    private class LockableModelListener implements ChangeListener, Serializable
    {
        public void stateChanged(ChangeEvent e)
        {
            if(!locked) fireStateChanged();
        }
    }

    public void repaint()
    {
        super.repaint();
    }

    public void repaint(Rectangle r)
    {
        super.repaint(r);
    }

	public void increment()
	{
		sliderModel.setValue(sliderModel.getValue()+1);
	}

	public void decrement()
	{
		sliderModel.setValue(sliderModel.getValue()-1);
	}
}
