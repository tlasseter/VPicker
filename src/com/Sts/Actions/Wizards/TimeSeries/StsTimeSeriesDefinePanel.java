package com.Sts.Actions.Wizards.TimeSeries;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsTimeSeriesDefinePanel extends JPanel implements ActionListener
{
    private StsTimeSeriesWizard wizard;
    private StsTimeSeriesDefine wizardStep;
    private StsModel model = null;

    JCheckBox cycleChkBox = new JCheckBox();
    JCheckBox cummChkBox = new JCheckBox();
  JPanel jPanel1 = new JPanel();
  JRadioButton actualRadio = new JRadioButton();
  JRadioButton artificialRadio = new JRadioButton();
  JTextField artificialTimeTxt = new JTextField();
  JLabel actualTimeLabel = new JLabel();
  JLabel jLabel1 = new JLabel();
  JLabel jLabel2 = new JLabel();
  JLabel jLabel3 = new JLabel();
  ButtonGroup actualArtificialGrp = new ButtonGroup();
  GridBagLayout gridBagLayout1 = new GridBagLayout();
  GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsTimeSeriesDefinePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsTimeSeriesWizard)wizard;
        this.wizardStep = (StsTimeSeriesDefine)wizardStep;
        this.model = wizard.getModel();

        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {

    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);
        this.setMinimumSize(new Dimension(200, 300));

        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout1);
        actualRadio.setSelected(true);
        actualRadio.setText("Actual Elapsed:");
        artificialRadio.setText("Artificial:");
//        artificialTimeTxt.setText("jTextField1");
        actualTimeLabel.setText("actualTime");
        actualRadio.addActionListener(this);
        artificialRadio.addActionListener(this);
        jLabel1.setText("minutes");
        jLabel2.setText("minutes");
        jLabel3.setFont(new java.awt.Font("Dialog", 1, 12));
        jLabel3.setHorizontalAlignment(SwingConstants.CENTER);
        jLabel3.setText("Playback Speed");
        jPanel1.add(jLabel3,  new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(1, 33, 0, 28), 58, 8));
        jPanel1.add(jLabel2,  new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 10, 7), 0, 5));
        jPanel1.add(actualRadio,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 4, 0, 0), 0, 0));
        jPanel1.add(artificialRadio,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 37, 10, 0), 0, -3));
        jPanel1.add(artificialTimeTxt,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 10, 8), -7, 2));
        jPanel1.add(jLabel1,  new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 7), 0, 5));
        jPanel1.add(actualTimeLabel,  new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 7, 5));
        actualArtificialGrp.add(actualRadio);
        actualArtificialGrp.add(artificialRadio);
    this.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 6, 4), 0, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == actualRadio)
        {
            artificialTimeTxt.setEnabled(false);
        }
        else if(source == artificialRadio)
        {
            artificialTimeTxt.setEnabled(true);
        }
    }

    public byte getTimeType()
    {
        if(actualRadio.isSelected())
            return wizard.ACTUAL_TIME;
        else
            return wizard.ARTIFICIAL_TIME;
    }

    public int getDurationMinutes()
    {
        if(getTimeType() == wizard.ACTUAL_TIME)
            return -1;
        else
            return new Integer(artificialTimeTxt.getText()).intValue();
    }
}
