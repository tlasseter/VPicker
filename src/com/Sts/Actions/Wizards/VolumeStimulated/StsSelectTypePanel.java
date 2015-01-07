package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.Well.StsFileType;
import com.Sts.Actions.Wizards.Well.StsWellWizard;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;

import javax.swing.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectTypePanel extends StsJPanel {
    private StsVolumeStimulatedWizard wizard;
    private StsSelectType wizardStep;

    JRadioButton fractureRadio = new JRadioButton();
    JRadioButton volumeRadio = new JRadioButton();
    ButtonGroup typeGroup = new ButtonGroup();
    StsGroupBox typeBox = new StsGroupBox("Select Computation Type");

    public StsSelectTypePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsVolumeStimulatedWizard)wizard;
        this.wizardStep = (StsSelectType)wizardStep;

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
        fractureRadio.setText("Use Fracture Sets");
        fractureRadio.setToolTipText("Use interpreted fracture sets to compute the volume stimulated.");
        volumeRadio.setText("Use Seismic Volume");
        volumeRadio.setToolTipText("Use the filled cells in a selected volume to compute the volume stimulated");
        typeGroup.add(fractureRadio);
        typeGroup.add(volumeRadio);

        typeBox.add(fractureRadio);
        typeBox.add(volumeRadio);
        this.gbc.fill = gbc.HORIZONTAL;
        add(typeBox);

        // Always default to fracture set type.
        fractureRadio.setSelected(true);
    }

    public byte getType()
    {
        if(fractureRadio.isSelected()) return wizard.FRACTURES;
        if(volumeRadio.isSelected()) return wizard.VOLUME;
        return (byte)-1;
    }
}