
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI.Beans;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;


/** A fieldBean used currently in the objectPanel to display an uneditable colorscale.
 *  The bean objectType is an StsColorscale, i.e., we get/set colorscales into this
 *  bean as the panel holding the bean changes its panelObject as when we select different
 *  instances of a class which displays this fieldbean.
 */

public class StsColorscaleFieldBean extends StsFieldBean implements ActionListener
{
    private StsColorscale colorscale = null;
    private StsColorscalePanel colorscalePanel = new StsColorscalePanel(this);

    public StsColorscaleFieldBean()
    {
		this.add(colorscalePanel);
    }

    public StsColorscaleFieldBean(Class c, String fieldName)
    {
        super();
        super.initialize(c, fieldName);
		this.add(colorscalePanel);
    }

    public Component[] getBeanComponents() { return new Component[] { this }; }
//    public StsColorscalePanel getPanel() { return this; }

    public void addListener(ActionListener l) { }

    public Object getValueObject() { return getColorscale(); }
    public StsColorscale getColorscale() { return colorscale; }

    public void doSetValueObject(Object colorscaleObject)
    {
        if(!(colorscaleObject instanceof StsColorscale)) return;
        colorscale = (StsColorscale)colorscaleObject;
        colorscalePanel.setColorscale(colorscale);
        colorscalePanel.repaint();
    }

    public void actionPerformed(ActionEvent e)
    {
        System.out.println("Colorscale actionEvent: " + e.getSource().toString());
    }

	public String toString() {return StsParameters.NONE_STRING; }
	public Object fromString(String string) { return null; }

    static public void main(String[] args)
    {
        try
        {
            StsModel model = new StsModel();
//            StsSpectrum spectrum = new StsSpectrum("test", colors);
            CFBTest test = new CFBTest();
            StsColorscaleFieldBean colorscaleBean = new StsColorscaleFieldBean();
            colorscaleBean.initialize(test, "colorscale");

            JFrame frame = new JFrame();

//            StsColorscalePanel panel = colorscaleBean.getPanel();
            frame.getContentPane().add(colorscaleBean);

            frame.pack();
            frame.setVisible(true);
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}

class CFBTest
{
    StsColorscale colorscale = null;

        StsColor[] stsColors =
        {
            new StsColor(Color.blue),
            new StsColor(Color.cyan),
            new StsColor(Color.green),
            new StsColor(Color.magenta),
            new StsColor(Color.orange),
            new StsColor(Color.pink),
            new StsColor(Color.red),
            new StsColor(Color.yellow)
        };

    CFBTest()
    {
        StsSpectrum spectrum = StsSpectrum.constructor("test", stsColors, 255);
        colorscale = new StsColorscale("Test", spectrum, -100.0f, 100.0f);
    }

    public StsColorscale getColorscale() { return colorscale; }
    public void setColorscale(StsColorscale colorscale) { this.colorscale = colorscale; }
}

/*
    public void setValue(StsColorscale cscale)
    {
    	if( cscale != null )
        {
        	if( cscale.getSpectrum() != null )
	        	panel.setColors(cscale.getSpectrum().getColors());

            if( cscale.range() != null )
            	panel.setRange(cscale.range());

            panel.repaint();
        }
    }
*/
/*
    public void setValue()
    {
    	try
        {
            if(objectPanel == null) return;
            Object obj = objectPanel.getObject();
	        StsColorscale cscale = (StsColorscale) get.invoke(obj, new Object[0]);
    		setValue(cscale);
        }
        catch(Exception e) { e.printStackTrace(); }
    }
*/
