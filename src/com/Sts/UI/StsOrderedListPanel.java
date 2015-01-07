
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;


import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class StsOrderedListPanel extends JPanel implements ActionListener, ListSelectionListener
{
	GridBagLayout gridBagLayout = new GridBagLayout();
	JLabel titleLabel = new JLabel();
	JList list = new JList();
	JButton upButton = new JButton();
	JButton downButton = new JButton();
    Object[] items = null;

	public StsOrderedListPanel()
    {
		try { jbInit(); }
        catch(Exception e) { e.printStackTrace(); }
    }


    private void jbInit()
    {
		titleLabel.setText("Title:");
		upButton.setText("Up");
		downButton.setText("Down");
        upButton.setEnabled(false);
        downButton.setEnabled(false);
		upButton.addActionListener(this);
		downButton.addActionListener(this);
        list.addListSelectionListener(this);

        JScrollPane scrollPane = new JScrollPane();
        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
        scrollPane.getViewport().add(list);

        this.setLayout(gridBagLayout);
		this.add(titleLabel, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		this.add(scrollPane, new GridBagConstraints(0, 1, 1, 2, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
		this.add(upButton, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.SOUTH, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));
		this.add(downButton, new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(0, 5, 5, 5), 0, 0));
    }

    public void setTitle(String title) { titleLabel.setText(title); }
    public String getTitle() { return titleLabel.getText(); }

    public void setListRenderer(ListCellRenderer renderer)
    {
    	list.setCellRenderer(renderer);
    }

    public void setItems(Object[] items)
    {
    	this.items = items;
		list.setListData(items);
    }

    public Object[] getItems() { return items; }

    public Object[] getSelectedItems() { return list.getSelectedValues(); }
    public int[] getSelectedIndices() { return list.getSelectedIndices(); }
    public void setSelectedIndices(int[] selected) { list.setSelectedIndices(selected); }

    public void actionPerformed(ActionEvent e)
    {
    	if( e.getSource() == upButton ) moveUp();
    	else if( e.getSource() == downButton ) moveDown();

    }

    public void valueChanged(ListSelectionEvent e)
    {
    	int selectedIndex = list.getSelectedIndex();
    	if( list.getSelectedIndex() == 0 ) upButton.setEnabled(false);
        else upButton.setEnabled(true);
    	if( list.getSelectedIndex() == items.length-1 ) downButton.setEnabled(false);
        else downButton.setEnabled(true);
    }
    private void moveUp()
    {
        int index = list.getSelectedIndex();
        if (index==0) return;
        StsMath.swap(items, index-1, index);
        setItems(items);
        list.setSelectedIndex(index-1);
    }

    private void moveDown()
    {
        int index = list.getSelectedIndex();
        if( index == items.length-1)  return;
        StsMath.swap(items, index, index+1);
        setItems(items);
        list.setSelectedIndex(index+1);
    }



}
