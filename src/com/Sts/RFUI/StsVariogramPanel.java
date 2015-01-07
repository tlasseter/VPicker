

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;
public class StsVariogramPanel extends JPanel
{

    StsGroupBox groupBox1 = new StsGroupBox("Variogram Parameters");
    GridLayout gridLayout1 = new GridLayout();
    GridLayout gridLayout2 = new GridLayout();
    JLabel jLabel1 = new JLabel();
    JTextField maxRangeField = new JTextField();
    JLabel jLabel2 = new JLabel();
    JTextField minRangeField = new JTextField();
    JLabel jLabel3 = new JLabel();
    JTextField vertRangeField = new JTextField();
    JLabel jLabel4 = new JLabel();
    JTextField aziField = new JTextField();
    JLabel jLabel5 = new JLabel();
    JTextField sillField = new JTextField();
    JLabel jLabel6 = new JLabel();
    JTextField nuggetField = new JTextField();
    JLabel jLabel7 = new JLabel();
    JComboBox jComboBox1 = new JComboBox();


    String[] modelTypes = new String[]
    {
        "fBm",
        "fGn",
        "Linear",
        "Exponential",
        "Spherical",
        "Gaussian"
    };

    public StsVariogramPanel()
    {
        try
        {
            jbInit();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
        this.setLayout(gridLayout1);
        gridLayout2.setColumns(2);
        gridLayout2.setHgap(5);
        gridLayout2.setRows(0);
        gridLayout2.setVgap(2);
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Maximum Range");
        jLabel2.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel2.setText("Minimum Range");
        jLabel3.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel3.setText("Vertical Range");
        jLabel4.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel4.setText("Areal Azimuth");
        jLabel5.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel5.setText("Sill");
        jLabel6.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel6.setText("Nugget");
        jLabel7.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel7.setText("Model Type");
        for( int i=0; i<modelTypes.length; i++ ) jComboBox1.addItem(modelTypes[i]);

        this.add(groupBox1, null);
        groupBox1.add(jLabel1, null);
        groupBox1.add(maxRangeField, null);
        groupBox1.add(jLabel2, null);
        groupBox1.add(minRangeField, null);
        groupBox1.add(jLabel3, null);
        groupBox1.add(vertRangeField, null);
        groupBox1.add(jLabel4, null);
        groupBox1.add(aziField, null);
        groupBox1.add(jLabel5, null);
        groupBox1.add(sillField, null);
        groupBox1.add(jLabel6, null);
        groupBox1.add(nuggetField, null);
        groupBox1.add(jLabel7, null);
        groupBox1.add(jComboBox1, null);
    }


    public void setFloat(JTextField field, float f)
    {
        field.setText(Float.toString(f));
    }

    public Float getFloat(JTextField f)
    {
        return Float.valueOf(f.getText());
    }

    public void setMaxRange(float f) { setFloat(maxRangeField, f); }
    public Float getMaxRange() { return getFloat(maxRangeField); }
    public void setMinRange(float f) { setFloat(minRangeField, f); }
    public Float getMinRange() { return getFloat(minRangeField); }
    public void setVertRange(float f) { setFloat(vertRangeField, f); }
    public Float getVertRange() { return getFloat(vertRangeField); }
    public void setAzi(float f) { setFloat(aziField, f); }
    public Float getAzi() { return getFloat(aziField); }
    public void setSill(float f) { setFloat(sillField, f); }
    public Float getSill() { return getFloat(sillField); }
    public void setNugget(float f) { setFloat(nuggetField, f); }
    public Float getNugget() { return getFloat(nuggetField); }
    public void setModelType(int index) { jComboBox1.setSelectedIndex(index); }
    public int getModelType() { return jComboBox1.getSelectedIndex(); }


    public static void main(String[] args)
    {
        StsVariogramPanel StsVariogramPanel = new StsVariogramPanel();
        JFrame frame = new JFrame();
        frame.getContentPane().add(StsVariogramPanel);
        frame.pack();
        frame.setVisible(true);
    }
}
