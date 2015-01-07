

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsGroupCheckBox extends Panel implements ActionListener
{
    private EtchedRectangle box = new EtchedRectangle();
	JCheckBox checkBox = new JCheckBox();
	GridBagLayout gridBagLayout1 = new GridBagLayout();
	JPanel containerPanel = new JPanel();
	private transient Vector actionListeners = null;
    Font defaultFont = null;

	public StsGroupCheckBox()
	{
		try
		{
			jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	private void jbInit() throws Exception
	{
		setLayout(gridBagLayout1);
		checkBox.setText("Label");
		add(checkBox, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(0, 8, 0, 0), 0, 0));
		add(containerPanel, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));

        checkBox.addActionListener(this);
		validate();
	}

    public void paint(Graphics g)
    {
//    	System.out.println("Painting box...");
		box.paint(g);
        Dimension d = checkBox.getSize();
        g.setColor(getBackground());
        g.fillRect(5, 0, d.width+5, d.height);
		super.paint(g);
    }

    public void setBounds(int x, int y, int w, int h)
    {
		super.setBounds(x,y,w,h);
        FontMetrics fm   = checkBox.getFontMetrics(
                                checkBox.getFont());
        int         top  = getInsets().top + fm.getAscent();
        Dimension   size = getSize();

        box.setBounds(0, top, size.width-1, size.height-top-1);
	}

    public void setLabel(String title)
    {
    	checkBox.setText(title);
    }
	public String getLabel() { return checkBox.getText(); }

    public void setState(boolean selected, boolean grayed)
    {
    	setSelected(selected);
        setGrayed(grayed);
    }

    public void setGrayed(boolean b)
    {
        if( defaultFont == null )
            defaultFont = checkBox.getFont();

    	if( b )
	    	checkBox.setFont(new Font("Dialog", Font.ITALIC, defaultFont.getSize()));
        else
	    	checkBox.setFont(defaultFont);
    }

    public void setEnabled(boolean b)
    {
    	checkBox.setEnabled(b);
    }
    public void setSelected(boolean b)
    {
    	checkBox.setSelected(b);
        enableAll(b);
//        System.out.println("Checkbox set: " + b);
    }
    public boolean isSelected() { return checkBox.isSelected(); }

	public void setContainer(JPanel newContainer)
	{
		containerPanel = newContainer;
	}
    public JPanel getContainer() { return containerPanel; }

    public void actionPerformed(ActionEvent e)
    {
//        System.out.println("Checkbox action event: " + checkBox.isSelected());
		enableAll(checkBox.isSelected());
		setGrayed(false);
        fireActionPerformed(e);
    }

    public void enableAll(boolean b)
    {
    	Component[] components = containerPanel.getComponents();
        if( components != null )
        	for( int i=0; i<components.length; i++ )
		    	components[i].setEnabled(b);
	}

	protected void fireActionPerformed(ActionEvent e)
	{
		if (actionListeners != null)
		{
			Vector listeners = actionListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			  ((ActionListener) listeners.elementAt(i)).actionPerformed(e);
		}
	}

	public synchronized void removeActionListener(ActionListener l)
	{
		if (actionListeners != null && actionListeners.contains(l))
		{
			Vector v = (Vector) actionListeners.clone();
			v.removeElement(l);
			actionListeners = v;
		}
	}

	public synchronized void addActionListener(ActionListener l)
	{
		Vector v = actionListeners == null ? new Vector(2) : (Vector) actionListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			actionListeners = v;
		}
	}

    public void setActionCommand(String s) { checkBox.setActionCommand(s); }
    public String getActionCommand() { return checkBox.getActionCommand(); }


}

