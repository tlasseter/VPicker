package com.Sts.Actions.Wizards.AnalogueCube;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */
public class StsDefineTypePanel extends JPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;
    JPanel jPanel1 = new JPanel();
//    JCheckBox persistChk = new JCheckBox();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    StsGroupBox typeGroupBox = new StsGroupBox();
    JRadioButton realBtn = new JRadioButton("Real");
    JRadioButton complexBtn = new JRadioButton("Complex");

    ButtonGroup group = new ButtonGroup();

    public StsDefineTypePanel(StsWizard wizard, StsWizardStep wizardStep)
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

    }

    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout1);
        jPanel1.setLayout(gridBagLayout2);
//        persistChk.setText("Save Intermediate Volumes");
//        persistChk.setEnabled(true);

        typeGroupBox.setLabel("Analysis Type");
        typeGroupBox.setForeground(Color.gray);
        typeGroupBox.setFont(new java.awt.Font("Dialog", 0, 11));

        realBtn.setSelected(true);

        this.add(jPanel1,       new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 4, 0, 4), 0, 0));

        jPanel1.add(typeGroupBox,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 2, 2, 5), 3, 0));
        typeGroupBox.add(realBtn,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 2, 2, 5), 3, 0));
        typeGroupBox.add(complexBtn,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(5, 2, 5, 5), 3, 0));

//        jPanel1.add(persistChk,  new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
//            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(20, 2, 5, 5), 3, 0));

        group.add(realBtn);
        group.add(complexBtn);
    }

    public byte getType()
    {
        if(realBtn.isSelected())
            return StsAnalogueCubeWizard.REAL;
        else
            return StsAnalogueCubeWizard.COMPLEX;
    }
/*
    public boolean persistVolumes()
    {
        return persistChk.isSelected();
    }
*/
}
