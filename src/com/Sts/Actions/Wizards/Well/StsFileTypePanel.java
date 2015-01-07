package com.Sts.Actions.Wizards.Well;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;

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

public class StsFileTypePanel extends StsJPanel {
    private StsWellWizard wizard;
    private StsFileType wizardStep;

    JPanel jPanel1 = new JPanel();
    JRadioButton s2sRadio = new JRadioButton();
    JRadioButton geographixRadio = new JRadioButton();
    ButtonGroup typeGroup = new ButtonGroup();
    StsGroupBox typeBox = new StsGroupBox("Select Format Type");

    public StsFileTypePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsWellWizard)wizard;
        this.wizardStep = (StsFileType)wizardStep;

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


        s2sRadio.setText("S2S Formatted Wells (well-dev, well-ref...)");
        s2sRadio.setToolTipText("Wells to be loaded are in an S2S documented format");
        geographixRadio.setText("Geographix Formatted Wells (.wls)");
        geographixRadio.setToolTipText("Wells to be loaded are in a Geographix format");
        typeGroup.add(s2sRadio);
        typeGroup.add(geographixRadio);

        typeBox.add(s2sRadio);
        typeBox.add(geographixRadio);
        this.gbc.fill = gbc.HORIZONTAL;
        add(typeBox);

        // Always default to single well.
        s2sRadio.setSelected(true);
    }

    public byte getType()
    {
        if(s2sRadio.isSelected()) return wizard.S2S_WELLS;
        if(geographixRadio.isSelected()) return wizard.GEOGRAPHIX_WELLS;
        return (byte)-1;
    }
}
