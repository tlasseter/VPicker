package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
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

public class StsDefineWellSubVolumePanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    StsWell[] availableWells = null;
    StsWell well = null;
    JTextField cvName = new JTextField();
    JLabel cvLabel = new JLabel("Name:");
    private StsModel model = null;
    JPanel jPanel1 = new JPanel();
    JComboBox wellCombo = new JComboBox();
    JLabel jLabel1 = new JLabel();
    JLabel jLabel2 = new JLabel();
    JLabel jLabel3 = new JLabel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();
    JLabel jLabel4 = new JLabel();
    JTextField topZText = new JTextField();
    JTextField btmZText = new JTextField();
    JLabel jLabel5 = new JLabel();
    JLabel jLabel6 = new JLabel();
    JTextField radiusText = new JTextField();
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public StsDefineWellSubVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
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
        model = wizard.getModel();
        availableWells = (StsWell[])model.getCastObjectList(StsWell.class);
        for(int i=0; i<availableWells.length; i++)
        {
            wellCombo.addItem(availableWells[i].getName());
        }
    }

    void jbInit() throws Exception
    {
        jPanel1.setBorder(BorderFactory.createEtchedBorder());
        jPanel1.setLayout(gridBagLayout1);
        this.setLayout(gridBagLayout2);
        jLabel1.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel1.setText("Wells:");
        wellCombo.addActionListener(this);
        jLabel4.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel4.setText("Top Z:");
        topZText.setText("0.0");
        btmZText.setText("0.0");
        jLabel5.setText("Bottom Z:");
        jLabel5.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel6.setHorizontalAlignment(SwingConstants.RIGHT);
        jLabel6.setText("Radius:");
        radiusText.setText("0.0");
        cvName.setText("SubVolumeName");
        this.add(jPanel1,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(5, 5, 5, 5), 0, 0));

        jPanel1.add(cvLabel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 5, 0, 5), 0, 0));
        jPanel1.add(cvName,  new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 5), 0, 0));

        jPanel1.add(wellCombo,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 5), 0, 0));
        jPanel1.add(jLabel1,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));

        jPanel1.add(jLabel4,  new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));
        jPanel1.add(topZText,  new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 5), 0, 0));
        jPanel1.add(btmZText,  new GridBagConstraints(1, 3, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 5), 0, 0));
        jPanel1.add(jLabel5,  new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));
        jPanel1.add(radiusText,  new GridBagConstraints(1, 4, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(5, 5, 0, 5), 0, 0));
        jPanel1.add(jLabel6,  new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0
            ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 5, 0, 0), 0, 0));
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
        if(source == wellCombo)
        {
            well= availableWells[wellCombo.getSelectedIndex()];
        }
    }

    public StsWell getWell()
    {
        return well;
    }
    public float getTopZ()
    {
        return (new Float(topZText.getText()).floatValue());
    }

    public float getBottomZ()
    {
        return (new Float(btmZText.getText()).floatValue());
    }

    public float getRadius()
    {
        return (new Float(radiusText.getText()).floatValue());
    }

    public String getSubVolumeName()
    {
        return cvName.getText();
    }

}
