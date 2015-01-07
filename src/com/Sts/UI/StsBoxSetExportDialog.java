package com.Sts.UI;

import com.Sts.DBTypes.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsBoxSetExportDialog extends JDialog implements ActionListener
{
    JPanel jPanel1 = new JPanel();
    ButtonGroup resGrp = new ButtonGroup();

    public final static byte XML = 0;
    public final static byte ASCII = 1;

    public final static byte CENTERS = 0;
    public final static byte BOXES = 1;

    public final static byte PROCESS = 0;
    public final static byte CANCELED = 1;

    byte mode = PROCESS;
    byte format = ASCII;
    byte scope = CENTERS;
    boolean outAmplitude = false;

    StsBoxSetSubVolume subVolume = null;

    JButton processBtn = new JButton();
    JButton cancelBtn = new JButton();
    JPanel jPanel2 = new JPanel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    TitledBorder titledBorder1;
    JPanel jPanel3 = new JPanel();
    JRadioButton asciiRadio = new JRadioButton();
    JRadioButton centersRadio = new JRadioButton();
    JLabel fileViewLbl = new JLabel();
    JRadioButton xmlRadio = new JRadioButton();
    JRadioButton boxesRadio = new JRadioButton();
    Frame frame = null;
    ButtonGroup group = new ButtonGroup();
    JPanel jPanel4 = new JPanel();
    JLabel jLabel1 = new JLabel();
    JPanel jPanel5 = new JPanel();
    JLabel jLabel2 = new JLabel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    GridBagLayout gridBagLayout3 = new GridBagLayout();
  JCheckBox amplitudeChkBox = new JCheckBox();
  GridBagLayout gridBagLayout4 = new GridBagLayout();
  GridBagLayout gridBagLayout5 = new GridBagLayout();
  GridBagLayout gridBagLayout6 = new GridBagLayout();

    public StsBoxSetExportDialog(Frame frame, String title, boolean modal, StsBoxSetSubVolume subVolume)
    {
        super(frame, title, modal);
        this.frame = frame;
        this.subVolume = subVolume;
        try
        {
            jbInit();
            pack();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public StsBoxSetExportDialog()
    {
        this(null, "", false, null);
    }

    private void jbInit() throws Exception
    {
        titledBorder1 = new TitledBorder("");
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout5);
        this.getContentPane().setLayout(gridBagLayout6);
        processBtn.setText("Process");
        cancelBtn.setText("Cancel");
        jPanel2.setBorder(BorderFactory.createEtchedBorder());
        jPanel2.setLayout(gridBagLayout1);
        processBtn.addActionListener(this);
        cancelBtn.addActionListener(this);
        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setLayout(gridBagLayout4);
        asciiRadio.addActionListener(this);
        asciiRadio.setText("ASCII");
        centersRadio.setSelected(true);
    centersRadio.setText("Centers Only");
        centersRadio.addActionListener(this);
        fileViewLbl.setFont(new java.awt.Font("Serif", 1, 14));
        fileViewLbl.setHorizontalAlignment(SwingConstants.CENTER);
        fileViewLbl.setText("Box Set Sub-PostStack3d Export");
        asciiRadio.setSelected(true);
        xmlRadio.setText("XML");
        xmlRadio.addActionListener(this);
        boxesRadio.setText("Complete Boxes");
        boxesRadio.addActionListener(this);
        jLabel1.setForeground(Color.black);
        jPanel4.setLayout(gridBagLayout3);
        jLabel1.setFont(new java.awt.Font("Dialog", 1, 11));
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Scope:");
        jLabel2.setFont(new java.awt.Font("Dialog", 1, 11));
        jLabel2.setText("Format:");
        jPanel5.setLayout(gridBagLayout2);
        amplitudeChkBox.setText("Include PostStack3d Info");
        amplitudeChkBox.addActionListener(this);
        amplitudeChkBox.setSelected(false);
    jPanel2.add(processBtn,   new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 1, 3, 0), 0, 0));
        jPanel2.add(cancelBtn,    new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(2, 46, 2, 2), 4, 0));
        jPanel1.add(fileViewLbl,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 1, 0, 0), 154, 4));
        jPanel1.add(jPanel3,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 1, 0, 0), 0, 1));
        jPanel4.add(jLabel1,   new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, -3, 2, 0), 18, 10));
        jPanel4.add(centersRadio,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 23, 2, 0), 14, 2));
        jPanel4.add(boxesRadio,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 17, 2, 18), 0, 2));
    jPanel3.add(amplitudeChkBox, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 198, 0, 2), 9, 0));
    jPanel1.add(jPanel2, new GridBagConstraints(0, 2, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 1, 0, 0), 133, 3));
    jPanel3.add(jPanel5,  new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, -1, 0, 2), 0, -5));
        jPanel5.add(jLabel2,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 17, 7, 0), 0, 10));
        jPanel5.add(xmlRadio,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 22, 7, 0), 0, 2));
        jPanel5.add(asciiRadio,  new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(1, 73, 7, 72), 0, 2));
    jPanel3.add(jPanel4, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, -1, 0, 2), 6, 0));
        resGrp.add(xmlRadio);
        resGrp.add(asciiRadio);
    group.add(centersRadio);
    group.add(boxesRadio);
    this.getContentPane().add(jPanel1, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.SOUTH, GridBagConstraints.BOTH, new Insets(0, 3, 5, 5), 0, 7));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == xmlRadio)
        {
            format = XML;
            amplitudeChkBox.setEnabled(false);
        }
        else if(source == centersRadio)
        {
            scope = CENTERS;
        }
        else if(source == boxesRadio)
        {
            scope = BOXES;
        }
        else if(source == asciiRadio)
        {
            format = ASCII;
            amplitudeChkBox.setEnabled(true);
        }
        else if(source == processBtn)
        {
            mode = PROCESS;
            setVisible(false);
        }
        else if(source == cancelBtn)
        {
            mode = CANCELED;
            setVisible(false);
        }
        else if(source == amplitudeChkBox)
        {
            if (amplitudeChkBox.isSelected())
                outAmplitude = true;
            else
                outAmplitude = false;
        }
    }

    public byte getMode() { return mode;  }
    public float getFormat()
    {
        return format;
    }
    public byte getScope()
    {
        return scope;
    }

    public boolean getAmplitude()
    {
        return outAmplitude;
    }
}
