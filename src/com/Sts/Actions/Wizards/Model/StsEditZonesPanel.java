
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Model;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;

public class StsEditZonesPanel extends JPanel implements ListSelectionListener, ActionListener, FocusListener
{
    static final Object[] columns = { "Name ", "Type ", "Layer " };

	GridBagLayout gridBagLayout = new GridBagLayout();
	Border border = BorderFactory.createEtchedBorder();

	JLabel titleLabel = new JLabel();
	DefaultTableModel tableModel = new DefaultTableModel();
	JTable table = new JTable(tableModel);
	JScrollPane scrollPane = new JScrollPane(table);

    JPanel fieldPanel = new JPanel();
    JLabel nameLabel = new JLabel();
    JLabel typeLabel = new JLabel();
    JLabel nzLabel = new JLabel();
	JTextField nameField = new JTextField();
	JTextField nzField = new JTextField();
	JComboBox typeField = new JComboBox();

    StsZone[] zones = null;
    int[] selected = null;

	public StsEditZonesPanel()
    {
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
		typeLabel.setText("Type");
		nzLabel.setText("Layer");
        for( int i=0; i<StsZone.subZoneTypeStrings.length; i++ )
			typeField.addItem(StsZone.subZoneTypeStrings[i]);


        nameField.addActionListener(this);
        typeField.addActionListener(this);
        nzField.addActionListener(this);

        nameField.addFocusListener(this);
        typeField.addFocusListener(this);
        nzField.addFocusListener(this);

        addColumns(columns);
        table.getSelectionModel().addListSelectionListener(this);
		fieldPanel.add(nameLabel, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.SOUTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
		fieldPanel.add(typeLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
		fieldPanel.add(nzLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
		fieldPanel.add(nameField, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
		fieldPanel.add(nzField, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
		fieldPanel.add(typeField, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));

		this.setLayout(gridBagLayout);
		this.add(titleLabel, new GridBagConstraints(0, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 0, 0));
		this.add(scrollPane, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 5, 5, 5), 0, 0));
		this.add(fieldPanel, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 5, 5), 0, 0));
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

    public void addRows(Object[][] rows)
    {
    	int nRows = rows == null ? 0 : rows.length;
        for( int i=0; i<nRows; i++ )
			addRow(rows[i]);
//        table.sizeColumnsToFit(false);
    }

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
		if( zones == null ) return;
        if( selected != null && selected.length > 0 )
        {
        	int index = selected[0];
            if( index >= 0 && index < zones.length )
            {
				nameField.setText(zones[index].getName());
				typeField.setSelectedItem(zones[index].getSubZoneTypeStr());
				nzField.setText(Integer.toString(zones[index].getNSubZones()));
            }
        }
    }

    private void refreshList()
    {
		int[] indices = getSelectedIndices();
    	removeAllRows();
        int nZones = zones == null ? 0 : zones.length;
		Object[] row = new Object[3];
        for( int i=0; i<nZones; i++ )
        {
			row[0] = new String(zones[i].getName());
            row[1] = new String(zones[i].getSubZoneTypeStr());
            row[2] = new String(Integer.toString(zones[i].getNSubZones()));
            addRow(row);
        }
		setSelectedIndices(indices);
    }

    public void setZones(StsZone[] zones)
    {
    	this.zones = zones;
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
	    	zones[selected[0]].setName(text);
    }

    private void changeType(String text)
    {
    	for( int i=0; i<selected.length; i++ )
	    	zones[selected[i]].setSubZoneTypeStr(text);
    }

    private void changeNz(int nz)
    {
    	boolean beeped = false;
    	for( int i=0; i<selected.length; i++ )
        {
	    	zones[selected[i]].setNSubZones(nz);
	    	if( !beeped && zones[selected[i]].getNSubZones() != nz )
            {
				StsToolkit.beep();
                beeped = true;
            }
        }
    }

    private void updateField(Object field)
    {
		if( zones == null || selected == null ) return;
    	if( field == nameField )
        	changeName(nameField.getText());
    	else if( field == typeField )
        	changeType((String)typeField.getSelectedItem());
    	else if( field == nzField )
        {
        	String text = nzField.getText();
            text.trim();
            int nz = text == null || text.length() == 0 ? 0 : Integer.parseInt(text);
        	changeNz(nz);
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
