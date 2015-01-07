

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
import java.awt.event.*;
import java.io.*;

public class StsSimulationFileChoicePanel extends JPanel implements ActionListener
{
    GridBagLayout gridBagLayout = new GridBagLayout();
    GridLayout gridLayout1 = new GridLayout();

    JLabel jLabel1 = new JLabel();
    JTextField jTextField1 = new JTextField();
    JButton jButton1 = new JButton();

	ButtonGroup group = new ButtonGroup();
    StsGroupBox groupBox1 = new StsGroupBox("Format");
    JRadioButton jRadioButton1 = new JRadioButton();
    JRadioButton jRadioButton2 = new JRadioButton();

    File currentDir = new File(".");

    public StsSimulationFileChoicePanel()
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
        jLabel1.setText("Filename");
        jTextField1.setText("default");
        jButton1.setText("Browse...");
        jButton1.addActionListener(this);

        gridLayout1.setColumns(1);
        gridLayout1.setRows(0);
        jRadioButton1.setSelected(true);
        jRadioButton1.setText("Eclipse");
        jRadioButton2.setText("RESCUE");
        group.add(jRadioButton1);
        group.add(jRadioButton2);

        this.add(jLabel1, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 5, 0, 5), 0, 0));
        this.add(jTextField1, new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(10, 5, 0, 5), 0, 0));
        this.add(jButton1, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 0, 0, 0), 0, 0));
        this.add(groupBox1, new GridBagConstraints(0, 1, 2, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(50, 0, 0, 0), 0, 0));
        groupBox1.add(jRadioButton1, null);
        groupBox1.add(jRadioButton2, null);
    }

    public void setCurrentDir(File dir) { currentDir = dir; }

    public void actionPerformed(ActionEvent e)
    {
        JFileChooser chooser = new JFileChooser(currentDir);
        chooser.showDialog(this, "Select");
        if( chooser.getSelectedFile() != null )
        {
            File file = chooser.getSelectedFile();
            jTextField1.setText(file.getPath());
        }

    }

    public int getFormat()
    {
        if( jRadioButton1.isSelected() ) return 1;
        if( jRadioButton2.isSelected() ) return 2;
        return 0;
    }
    public String getFilename()
    {
        return jTextField1.getText();
    }

    public static void main(String[] args)
    {
        StsSimulationFileChoicePanel panel = new StsSimulationFileChoicePanel();
        JFrame frame = new JFrame();
        frame.getContentPane().add(panel);
        panel.setPreferredSize(new Dimension(300,150));
        frame.pack();
        frame.setVisible(true);
    }

}
