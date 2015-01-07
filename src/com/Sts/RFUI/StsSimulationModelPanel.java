

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

public class StsSimulationModelPanel extends JPanel
{
    GridBagLayout gridBagLayout = new GridBagLayout();

    JLabel yLabel = new JLabel();
    JTextField kxField = new JTextField();

    StsGroupBox groupBox1 = new StsGroupBox("Multipliers");
    JLabel permLabel = new JLabel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JLabel xLabel = new JLabel();
    JTextField kyField = new JTextField();
    JTextField kzField = new JTextField();
    JLabel zLabel = new JLabel();
    JComboBox permModelList = new JComboBox();
    JLabel porosityLabel = new JLabel();
    JComboBox porosityModelList = new JComboBox();
    JLabel jLabel1 = new JLabel();

    public StsSimulationModelPanel()
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
        this.setLayout(gridBagLayout);

        porosityLabel.setText("Porosity");
        permLabel.setText("Permeability");
        xLabel.setText("X");
        yLabel.setText("Y");
        zLabel.setText("Z");

        jLabel1.setText("Select property models");
        kxField.setText("1.0");
        kxField.setHorizontalAlignment(SwingConstants.RIGHT);
        kyField.setText("1.0");
        kyField.setHorizontalAlignment(SwingConstants.RIGHT);
        kzField.setText("0.1");
        kzField.setHorizontalAlignment(SwingConstants.RIGHT);
        add(porosityLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 10, 10), 0, 0));
        add(porosityModelList, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 10, 0), 0, 0));
        add(permLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 10, 10), 0, 0));
        add(permModelList, new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 10, 0), 0, 0));

        add(groupBox1, new GridBagConstraints(0, 3, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        groupBox1.add(xLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 10), 0, 0));
        groupBox1.add(yLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 10), 0, 0));
        groupBox1.add(zLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 10), 0, 0));

        groupBox1.add(kyField, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        groupBox1.add(kxField, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        groupBox1.add(kzField, new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        this.add(jLabel1, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 30, 0), 0, 0));

    }

    public void setPropertyModels(String[] models)
    {
        int n = models == null ? 0 : models.length;
        if( porosityModelList.getItemCount() > 0 )
        {
            porosityModelList.removeAllItems();
            permModelList.removeAllItems();
        }
        for( int i=0; i<n; i++ )
        {
            porosityModelList.addItem(models[i]);
            permModelList.addItem(models[i]);
        }
    }

    public int getPorosityModelIndex()
    {
        return porosityModelList.getSelectedIndex();
    }
    public int getPermeabilityModelIndex()
    {
        return permModelList.getSelectedIndex();
    }

    public void setDouble(JTextField field, double d)
    {
        field.setText(Double.toString(d));
    }

    public double getDouble(JTextField f)
    {
        return Double.valueOf(f.getText()).doubleValue();
    }

    public void setKx(double d) { setDouble(kxField, d); }
    public void setKy(double d) { setDouble(kyField, d); }
    public void setKz(double d) { setDouble(kzField, d); }

    public double getKx() { return getDouble(kxField); }
    public double getKy() { return getDouble(kyField); }
    public double getKz() { return getDouble(kzField); }


    public static void main(String[] args)
    {
        StsSimulationModelPanel panel = new StsSimulationModelPanel();
        JFrame frame = new JFrame();
        frame.getContentPane().add(panel);
        panel.setPreferredSize(new Dimension(300,150));
        frame.pack();
        frame.setVisible(true);
    }

}
