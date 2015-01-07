
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Model;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;


public class StsSurfaceAutopickPanel extends JPanel implements ListSelectionListener, ActionListener, FocusListener
{
	StsModel model;
    static final Object[] columns = { "Name ", "Min Correl Filter " };

	GridBagLayout gridBagLayout = new GridBagLayout();
	Border border = BorderFactory.createEtchedBorder();

	JLabel titleLabel = new JLabel();
	DefaultTableModel tableModel = new DefaultTableModel();
	JTable table = new JTable(tableModel);
	JScrollPane scrollPane = new JScrollPane(table);

    JPanel fieldPanel = new JPanel();
    JLabel nameLabel = new JLabel();
    JLabel correlLabel = new JLabel();
	JTextField nameField = new JTextField();
	JTextField correlField = new JTextField();

    StsModelSurface[] surfaces = null;
	StsHorpick[] horpicks = null;
    int[] selected = null;

	public StsSurfaceAutopickPanel(StsModel model)
    {
		this.model = model;
		try { jbInit(); }
        catch(Exception e) { e.printStackTrace(); }
    }

    private void jbInit()
    {
    	table.setShowVerticalLines(false);
		table.setShowHorizontalLines(false);
		table.setRequestFocusEnabled(false);
		scrollPane.setBorder(border);
		fieldPanel.setLayout(gridBagLayout);
		fieldPanel.setBorder(border);
		titleLabel.setText("Title");
		nameLabel.setText("Name");
		correlLabel.setText("Min Correl Filter");


        nameField.setEnabled(false);
        correlField.addActionListener(this);

        correlField.addFocusListener(this);

		addColumns(columns);
		table.getSelectionModel().addListSelectionListener(this);

		this.setLayout(gridBagLayout);
		this.add(titleLabel,  new GridBagConstraints(0, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 0, 0));
		this.add(scrollPane,   new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.EAST, GridBagConstraints.VERTICAL, new Insets(5, 5, 5, 5), 200, 0));
		this.add(fieldPanel, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
			,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(5, 5, 5, 5), 0, 0));

		fieldPanel.add(nameLabel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
		fieldPanel.add(correlLabel,   new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
				,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		fieldPanel.add(nameField,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 90, 0));
		fieldPanel.add(correlField,   new GridBagConstraints(1, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.NORTH, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 90, 0));

    }

    public void setTitle(String title) { titleLabel.setText(title); }

    public void addColumns(Object[] headers)
    {
    	int nCols = headers == null ? 0 : headers.length;
        for( int i=0; i<nCols; i++ )
	      	tableModel.addColumn(headers[i]);
    }

    public void addColumn(Object header)
    {
      	tableModel.addColumn(header);
    }
/*
    public void addRows(Object[][] rows)
    {
    	int nRows = rows == null ? 0 : rows.length;
        for( int i=0; i<nRows; i++ )
			addRow(rows[i]);
//        table.sizeColumnsToFit(false);
    }
*/
    public void addRow(Object[] row)
    {
		tableModel.addRow(row);
    }

    public int[] getSelectedIndices() { return table.getSelectedRows(); }
    public void setSelectedIndices(int[] indices)
    {
    	DefaultListSelectionModel model = (DefaultListSelectionModel) table.getSelectionModel();
    	if( indices == null ) model.clearSelection();
        else
        {
        	for( int i=0; i<indices.length; i++ )
            	model.addSelectionInterval(indices[i], indices[i]);
        }
    }

    public void removeAllRows()
    {
		tableModel.setNumRows(0);
    }

    private void refreshFields()
    {
		if( surfaces == null ) return;
        if( selected != null && selected.length > 0 )
        {
        	int index = selected[0];
            if( index >= 0 && index < surfaces.length )
            {
				nameField.setText(surfaces[index].getName());
				StsHorpick horpick = horpicks[index];
				if(horpick == null)
				{
					correlField.setText("0.0");
					correlField.setEnabled(false);
				}
				else
				{
					float minCorrel = horpick.getMinCorrelFilter();
					correlField.setText(Float.toString(minCorrel));
					correlField.setEnabled(true);
				}
            }
        }
    }

	private float getAutopickMinCorrel(StsModelSurface modelSurface)
	{
		StsSurface surface = modelSurface.getOriginalSurface();
		StsHorpickClass horpickClass = (StsHorpickClass) model.getStsClass(StsHorpick.class);
		StsHorpick horpick = horpickClass.getHorpickWithSurface(surface);
		if(horpick == null) return 0.0f;
		return horpick.getMinCorrelFilter();
	}

	private void setAutopickMinCorrel(StsModelSurface modelSurface, float minCorrelFilter)
	{
		StsSurface surface = modelSurface.getOriginalSurface();
		StsHorpickClass horpickClass = (StsHorpickClass) model.getStsClass(StsHorpick.class);
		StsHorpick horpick = horpickClass.getHorpickWithSurface(surface);
		if(horpick == null) return;
		horpick.setMinCorrelFilter(minCorrelFilter);
	}

    private void refreshList()
    {
		int[] indices = getSelectedIndices();
    	removeAllRows();
        int nSurfaces = surfaces == null ? 0 : surfaces.length;
		Object[] row = new Object[2];
        for( int i=0; i<nSurfaces; i++ )
        {
			row[0] = new String(surfaces[i].getName());
			float minCorrel = getAutopickMinCorrel(surfaces[i]);
            row[1] = new String(Float.toString(minCorrel));
            addRow(row);
        }
		setSelectedIndices(indices);
    }

    public void setSurfaces(StsModelSurface[] surfaces)
    {
    	this.surfaces = surfaces;
		horpicks = new StsHorpick[surfaces.length];
		StsHorpickClass horpickClass = (StsHorpickClass) model.getStsClass(StsHorpick.class);
		for(int n = 0; n < surfaces.length; n++)
		{
			StsSurface originalSurface = surfaces[n].getOriginalSurface();
			horpicks[n] = horpickClass.getHorpickWithSurface(originalSurface);
		}
		refreshList();
    }

    public void valueChanged(ListSelectionEvent e)
    {
        selected = getSelectedIndices();
		refreshFields();
    }

    private void changeName(String text)
    {
    	if( selected.length > 0 )
	    	surfaces[selected[0]].setName(text);
    }

    private void changeCorrel(String text)
    {
    	for( int i=0; i<selected.length; i++ )
		{
			float minCorrelValue = Float.parseFloat(text);
			setAutopickMinCorrel(surfaces[selected[i]], minCorrelValue);
		}
    }

    private void updateField(Object field)
    {
		if( surfaces == null || selected == null ) return;
    	if( field == correlField )
		{
        	String text = correlField.getText();
            text.trim();
        	changeCorrel(text);
		}
        refreshList();
    }

	public void actionPerformed(ActionEvent e)
	{
    	updateField(e.getSource());
        refreshFields();
	}

	public void focusGained(FocusEvent e)
	{
	}

	public void focusLost(FocusEvent e)
	{
    	updateField(e.getSource());
        refreshFields();
	}
}
