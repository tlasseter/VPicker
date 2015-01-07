
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Model;

import com.Sts.DBTypes.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;

public class StsEditFaultGapsPanel extends JPanel
implements ListSelectionListener, ActionListener, FocusListener
{
    static final Object[] columns = { "Name ", "Left ", "Right " };

	GridBagLayout gridBagLayout = new GridBagLayout();
	Border border = BorderFactory.createEtchedBorder();

	JLabel titleLabel = new JLabel();
	DefaultTableModel tableModel = new DefaultTableModel();
	JTable table = new JTable(tableModel);
	JScrollPane scrollPane = new JScrollPane(table);

    JPanel fieldPanel = new JPanel();
    JLabel nameLabel = new JLabel();
    JLabel leftGapLabel = new JLabel();
    JLabel rightGapLabel = new JLabel();
	JTextField nameField = new JTextField();
	JTextField leftGapField = new JTextField();
	JTextField rightGapField = new JTextField();

    StsMainObject[] sections = null;
    int[] selected = null;

	public StsEditFaultGapsPanel()
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
		leftGapLabel.setText("Left");
		rightGapLabel.setText("Right");

        nameField.addActionListener(this);
        leftGapField.addActionListener(this);
        rightGapField.addActionListener(this);

        nameField.addFocusListener(this);
        leftGapField.addFocusListener(this);
        rightGapField.addFocusListener(this);

        addColumns(columns);
        table.getSelectionModel().addListSelectionListener(this);

		this.setLayout(gridBagLayout);
		this.add(titleLabel, new GridBagConstraints(0, 0, 3, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 0), 0, 0));
		this.add(scrollPane, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
			,GridBagConstraints.EAST, GridBagConstraints.VERTICAL, new Insets(0, 0, 5, 5), 231, 0));
		this.add(fieldPanel, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
			,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(1, 0, 4, 3), -8, 0));

        fieldPanel.add(nameLabel, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.SOUTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
        fieldPanel.add(leftGapLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 32, 0));
        fieldPanel.add(rightGapLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 5, 5), 0, 0));
        fieldPanel.add(nameField, new GridBagConstraints(1, 0, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 4, 5, 0), 89, 0));
        fieldPanel.add(rightGapField, new GridBagConstraints(1, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.NORTHWEST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 90, 0));
        fieldPanel.add(leftGapField, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 89, 0));
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
		if( sections == null ) return;
        if( selected != null && selected.length > 0 )
        {
        	int index = selected[0];
            if( index >= 0 && index < sections.length )
            {
                StsSection section = (StsSection)sections[index];
				nameField.setText(section.getName());
				leftGapField.setText(Float.toString(section.getLeftGap()));
				rightGapField.setText(Float.toString(section.getRightGap()));
            }
        }
    }

    private void refreshList()
    {
		int[] indices = getSelectedIndices();
    	removeAllRows();
        int nSections = (sections == null) ? 0 : sections.length;
		Object[] row = new Object[3];
        for( int i=0; i<nSections; i++ )
        {
            StsSection section = (StsSection)sections[i];
			row[0] = new String(section.getName());
            row[1] = new String(Float.toString(section.getLeftGap()));
            row[2] = new String(Float.toString(section.getRightGap()));
            addRow(row);
        }
		setSelectedIndices(indices);
    }

    public void setSections(StsMainObject[] sections)
    {
    	this.sections = sections;
        refreshList();
    }

    public void valueChanged(ListSelectionEvent e)
    {
        selected = getSelectedIndices();
		refreshFields();
    }

    private void changeLeftGap(float gap)
    {
    	for( int i=0; i<selected.length; i++ )
        {
            StsSection section = (StsSection)sections[selected[i]];
	    	section.setLeftGap(gap);
        }
    }

    private void changeRightGap(float gap)
    {
    	for (int i=0; i<selected.length; i++)
        {
            StsSection section = (StsSection)sections[selected[i]];
	    	section.setRightGap(gap);
        }
    }

    private void updateField(Object field)
    {
		if (sections == null || selected == null) return;
        if (field == leftGapField)
        {
        	String text = leftGapField.getText();
            text.trim();
            float gap = (text == null || text.length() == 0) ? 0.0f
                    : Float.parseFloat(text);
        	changeLeftGap(gap);
        }
        else if (field == rightGapField)
        {
        	String text = rightGapField.getText();
            text.trim();
            float gap = (text == null || text.length() == 0) ? 0.0f
                    : Float.parseFloat(text);
        	changeRightGap(gap);
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
