
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.util.*;

public class StsLabeledSlider extends JPanel implements ChangeListener
{
    protected GridBagLayout gridBag = new GridBagLayout();
	protected JLabel label = new JLabel();
	protected JSlider slider = new JSlider();
	protected JLabel value = new JLabel();
	private transient Vector changeListeners;

	protected float minimum = 0f;
	protected float maximum = 0f;

    public StsLabeledSlider(String name)
    {
        super();
        setTitle(new String(name));
    }

    public StsLabeledSlider()
    {
        try { jbInit(); }
        catch(Exception e)
        {
	       	System.out.println("Exception in : StsLabeledSlider()\n" + e);
        }
    }

    private void jbInit() throws Exception
    {
        setLayout(gridBag);
      	Border border = BorderFactory.createLoweredBevelBorder();
		setBorder(border);

		label.setText("Title");
		label.setHorizontalAlignment(0);
		slider.setOrientation(JSlider.VERTICAL);
		slider.setMaximum(1000);
		value.setText("0");
		value.setHorizontalAlignment(0);
		value.setToolTipText("slider value");

        add(label, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        add(value, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        add(slider, new GridBagConstraints(0, 2, 1, 1, 0.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));

    	slider.addChangeListener(this);
    }

    public void stateChanged(ChangeEvent e)
    {
        setValueText(getValue());
        if( changeListeners != null )
        {
        	Enumeration list = changeListeners.elements();
            while( list.hasMoreElements() )
            {
            	ChangeListener l = (ChangeListener) list.nextElement();
	        	l.stateChanged(new ChangeEvent(this));
            }
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
		Vector v = changeListeners == null ? new Vector(2) : (Vector) changeListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			changeListeners = v;
		}
	}

	public void setTitle(String title)
    {
        setName(title);
        label.setText(title);
    }

	public String getTitle()
	{
		return getName();
	}

    public void setMinMax(float min, float max)
    {
        this.minimum = min;
        this.maximum = max;
        this.setValue(max);
    }

    public void setModel(BoundedRangeModel model)
    {
        slider.setModel(model);
    }
    public BoundedRangeModel getModel()
    {
        return slider.getModel();
    }
    public void setMajorTickSpacing(int n)
    {
        slider.setMajorTickSpacing(n);
    }
    public void setMinorTickSpacing(int n)
    {
        slider.setMinorTickSpacing(n);
    }
    public void setPaintTicks(boolean b)
    {
        slider.setPaintTicks(b);
    }
    public void setPaintLabels(boolean b)
    {
        slider.setPaintLabels(b);
    }
    public void setOrientation(int direction)
    {
        slider.setOrientation(direction);
    }
    public void setInverted(boolean b)
    {
        slider.setInverted(b);
    }

    public float getValue()
    {
		return minimum + ((maximum - minimum) * slider.getValue() / 1000.f);
    }

    public void setValue(float val)
    {
        // set the slider
		int newValue = (int) ((val - minimum) * 1000.f / (maximum - minimum));
        slider.getModel().setValue(newValue);

        setValueText(val);
	}

    public void setValueText(float val)
    {
        // set the string
        String valStr = new String();
        valStr += val;

        int pos = valStr.indexOf(".");
        if( pos > 0 )
	        value.setText(valStr.substring(0,pos));
        else
	        value.setText(valStr);
    }

	public void setMinimum(float newMinimum)
	{
        float value = getValue();
		minimum = newMinimum;
        setValue(value);
	}

	public float getMinimum()
	{
		return minimum;
	}

	public void setMaximum(float newMaximum)
	{
        float value = getValue();
		maximum = newMaximum;
        setValue(value);
	}

	public float getMaximum()
	{
		return maximum;
	}

	public static void main(String[] args)
	{
		StsLabeledSlider slider = new StsLabeledSlider();
        JFrame frame = new JFrame();

        slider.setMaximum(10f);
        slider.setValue(3f);
        slider.setTitle("My Slider");

        frame.getContentPane().add(slider);
        frame.pack();
        frame.setVisible(true);
	}
}
