
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Horizons;

import com.Sts.Types.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

public class StsModelGridPanel extends JPanel
	implements ActionListener, FocusListener
{
	GridBagLayout gridBagLayout = new GridBagLayout();
	StsGridDefinition def;
	JLabel minLabel = new JLabel();
	JTextField xOriginField = new JTextField();
	JTextField yOriginField = new JTextField();
	JLabel maxLabel = new JLabel();
	JTextField xSizeField = new JTextField();
	JTextField ySizeField = new JTextField();
	JLabel xLabel = new JLabel();
	JLabel yLabel = new JLabel();
	JTextField xIncField = new JTextField();
	JLabel incLabel = new JLabel();
	JTextField yIncField = new JTextField();
	JTextField nxField = new JTextField();
	JLabel nLabel = new JLabel();
	JTextField nyField = new JTextField();
	JLabel titleLabel = new JLabel();
	JLabel dummy = new JLabel();

	public StsModelGridPanel()
    {
		try { jbInit(); }
        catch(Exception e) { e.printStackTrace(); }
    }


    private void jbInit()
    {
		minLabel.setText("Minimum");
		maxLabel.setText("Maximum");
		xLabel.setText("X");
		xLabel.setHorizontalAlignment(0);
		yLabel.setHorizontalAlignment(0);
		incLabel.setText("Increment");
		nLabel.setText("Number");
		yLabel.setText("Y");

        xOriginField.addActionListener(this);
        yOriginField.addActionListener(this);
        xSizeField.addActionListener(this);
        ySizeField.addActionListener(this);
        xIncField.addActionListener(this);
        yIncField.addActionListener(this);
        nxField.addActionListener(this);
        nyField.addActionListener(this);

        xOriginField.addFocusListener(this);
        yOriginField.addFocusListener(this);
        xSizeField.addFocusListener(this);
        ySizeField.addFocusListener(this);
        xIncField.addFocusListener(this);
        yIncField.addFocusListener(this);
        nxField.addFocusListener(this);
        nyField.addFocusListener(this);

		this.setLayout(gridBagLayout);
        xOriginField.setBackground(SystemColor.control);
        xOriginField.setEnabled(false);
        yOriginField.setBackground(SystemColor.control);
        yOriginField.setEnabled(false);
        xSizeField.setBackground(SystemColor.control);
        xSizeField.setEnabled(false);
        ySizeField.setBackground(SystemColor.control);
        ySizeField.setEnabled(false);
        xIncField.setBackground(SystemColor.control);
        xIncField.setEnabled(false);
        yIncField.setBackground(SystemColor.control);
        yIncField.setEnabled(false);
        nxField.setBackground(SystemColor.control);
        nxField.setEnabled(false);
        nyField.setBackground(SystemColor.control);
        nyField.setEnabled(false);
        this.add(titleLabel, new GridBagConstraints(0, 0, 3, 1, 1.0, 0.5
            ,GridBagConstraints.NORTH, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		this.add(minLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 4, 2, 2), 0, 0));
		this.add(xOriginField, new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.HORIZONTAL, new Insets(3, 3, 3, 3), 0, 0));
		this.add(yOriginField, new GridBagConstraints(2, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(3, 3, 3, 3), 0, 0));
		this.add(maxLabel, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 4, 2, 2), 0, 0));
		this.add(xSizeField, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 3, 3, 3), 0, 0));
		this.add(ySizeField, new GridBagConstraints(2, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 3, 3, 3), 0, 0));
		this.add(xLabel, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		this.add(yLabel, new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
		this.add(xIncField, new GridBagConstraints(1, 4, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(3, 3, 3, 3), 0, 0));
		this.add(incLabel, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 4, 2, 2), 0, 0));
		this.add(yIncField, new GridBagConstraints(2, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 3, 3, 3), 0, 0));
		this.add(nxField, new GridBagConstraints(1, 5, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 3, 3, 3), 0, 0));
		this.add(nLabel, new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 4, 2, 2), 0, 0));
		this.add(nyField, new GridBagConstraints(2, 5, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 3, 3, 3), 0, 0));
		this.add(dummy, new GridBagConstraints(0, 6, 3, 2, 1.0, 0.5
            ,GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(0, 0, 0, 0), 0, 0));
    }

    public void setTitle(String title) { titleLabel.setText(title); }
    public StsGridDefinition getGridDefinition() { return def; }
    public void setGridDefinition(StsGridDefinition def)
    {
    	this.def = def;
        refreshFields();
    }

    public void refreshFields()
    {
    	if( def != null )
        {
			xOriginField.setText(String.valueOf(def.getXOrigin()));
			yOriginField.setText(String.valueOf(def.getYOrigin()));
			xSizeField.setText(String.valueOf(def.getXSize()));
			ySizeField.setText(String.valueOf(def.getYSize()));
			xIncField.setText(String.valueOf(def.getXInc()));
			yIncField.setText(String.valueOf(def.getYInc()));
			nxField.setText(String.valueOf(def.getNCols()));
			nyField.setText(String.valueOf(def.getNRows()));
            validate();
            repaint();
        }
    }

    private void updateField(JTextField field)
    {
        String text = field.getText();

        if( field == xOriginField )
        {
            double value = Double.valueOf(text).doubleValue();
            double xOrigin = def.getXOrigin();
            def.setXMin((float)(value-xOrigin));
        }
        else if( field == yOriginField )
        {
            double value = Double.valueOf(text).doubleValue();
            double yOrigin = def.getYOrigin();
            def.setYMin((float)(value - yOrigin));
        }
        // need to implement this
        else if( field == xSizeField )
        {
            float value = Float.valueOf(text).floatValue();
//            def.setXSize(value);
        }
        // need to implement this
        else if( field == ySizeField )
        {
            float value = Float.valueOf(text).floatValue();
//            def.setYSize(value);
        }
        else if( field == xIncField )
        {
            float value = Float.valueOf(text).floatValue();
            def.setXInc(value);
        }
        else if( field == yIncField )
        {
            float value = Float.valueOf(text).floatValue();
            def.setYInc(value);
        }
        else if( field == nxField )
        {
            int value = Integer.parseInt(text);
            def.setNCols(value);
        }
        else if( field == nyField )
        {
            int value = Integer.parseInt(text);
            def.setNRows(value);
        }
		refreshFields();
    }
    public void actionPerformed(ActionEvent e)
    {
    	if( def == null ) return;
		if( e.getSource() instanceof JTextField )
        	updateField((JTextField) e.getSource());
    }

	public void focusLost(FocusEvent e)
    {
    	if( e.getSource() instanceof JTextField )
        	updateField((JTextField) e.getSource());
    }
	public void focusGained(FocusEvent e) { }
}
