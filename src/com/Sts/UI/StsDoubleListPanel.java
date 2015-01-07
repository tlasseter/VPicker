
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsDoubleListPanel extends JPanel
	implements ActionListener, ListSelectionListener{
	GridBagLayout gridBagLayout = new GridBagLayout();
	JLabel titleLabel = new JLabel();
	JLabel rightLabel = new JLabel();
	JLabel leftLabel = new JLabel();
    JScrollPane leftScrollPane = new JScrollPane();
    JScrollPane rightScrollPane = new JScrollPane();
	JList leftList = new JList();
    JList rightList =  new JList();
	JButton leftButton = new JButton();
	JButton rightButton = new JButton();
    Object[] leftItems = null;
    Object[] rightItems = null;
	private transient Vector actionListeners = null;
	private transient Vector listSelectionListeners = null;

	public StsDoubleListPanel()
    {
		try { jbInit(); }
        catch(Exception e) { e.printStackTrace(); }
    }


    private void jbInit()
    {
		leftScrollPane.getViewport().add(leftList);
        rightScrollPane.getViewport().add(rightList);
		leftButton.setText("<");
		rightButton.setText(">");
		titleLabel.setText("Title:");
		rightLabel.setText("right");
		leftLabel.setText("left");
        leftButton.setEnabled(false);
        rightButton.setEnabled(false);
		leftButton.addActionListener(this);
		rightButton.addActionListener(this);
        leftList.addListSelectionListener(this);
        rightList.addListSelectionListener(this);
		this.setLayout(gridBagLayout);
		this.add(titleLabel, new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
		this.add(leftLabel, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		this.add(leftScrollPane, new GridBagConstraints(0, 2, 1, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
		this.add(leftButton, new GridBagConstraints(1, 2, 1, 1, 0.0, 1.0
            ,GridBagConstraints.SOUTH, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 2, 5), 0, 0));
		this.add(rightButton, new GridBagConstraints(1, 3, 1, 1, 0.0, 1.0
            ,GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(2, 5, 0, 5), 0, 0));
		this.add(rightScrollPane, new GridBagConstraints(2, 2, 1, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
		this.add(rightLabel, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));

    }

    public void setTitle(String title) { titleLabel.setText(title); }
    public String getTitle() { return titleLabel.getText(); }

    public void setLeftTitle(String title) { leftLabel.setText(title); }
    public void setRightTitle(String title) { rightLabel.setText(title); }
    public JButton getLeftButton() { return leftButton; }
    public JButton getRightButton() { return rightButton; }

    public void setLeftListRenderer(ListCellRenderer renderer)
    {
    	leftList.setCellRenderer(renderer);
    }
    public void setLeftItems(Object[] items)
    {
    	this.leftItems = items;
		leftList.setListData(items);
        leftList.repaint();
    }
    public Object[] getLeftItems() { return leftItems; }
    public JList getLeftList() { return leftList; }
    public Object[] getSelectedLeftItems() { return leftList.getSelectedValues(); }
    public int[] getSelectedLeftIndices() { return leftList.getSelectedIndices(); }
    public void clearLeftSelections() { leftList.clearSelection(); }
    public void setSelectedLeftIndices(int[] selected)
    {
    	leftList.setSelectedIndices(selected);
    }

	public void enableLeftButton(boolean b) { leftButton.setEnabled(b); }

    public Object[] getRightItems() { return rightItems; }
    public JList getRightList() { return rightList; }
    public Object[] getSelectedRightItems() { return rightList.getSelectedValues(); }
    public void clearRightSelections() { rightList.clearSelection(); }
    public void setSelectedRightIndices(int[] selected)
    {
        rightList.setSelectedIndices(selected);
    }
	public void enableRightButton(boolean b) { rightButton.setEnabled(b); }
    public void setRightListRenderer(ListCellRenderer renderer)
    {
    	rightList.setCellRenderer(renderer);
    }
    public void setRightItems(Object[] items)
    {
    	this.rightItems = items;
		rightList.setListData(items);
        rightList.repaint();
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

	public void actionPerformed(ActionEvent e)
	{
    	fireActionPerformed(e);
	}

	public synchronized void removeListSelectionListener(ListSelectionListener l)
	{
		if (listSelectionListeners != null && listSelectionListeners.contains(l))
		{
			Vector v = (Vector) listSelectionListeners.clone();
			v.removeElement(l);
			listSelectionListeners = v;
		}
	}

	public synchronized void addListSelectionListener(ListSelectionListener l)
	{
		Vector v = listSelectionListeners == null ? new Vector(2) : (Vector) listSelectionListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			listSelectionListeners = v;
		}
	}

	protected void fireValueChanged(ListSelectionEvent e)
	{
		if (listSelectionListeners != null)
		{
			Vector listeners = listSelectionListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			  ((ListSelectionListener) listeners.elementAt(i)).valueChanged(e);
		}
	}

	public void valueChanged(ListSelectionEvent e)
	{
    	fireValueChanged(e);
	}
}
